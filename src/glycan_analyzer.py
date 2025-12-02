"""
Glycan Analyzer Module
======================

This module calculates glycan (carbohydrate) requirements for HMPV glycoproteins
for inclusion in the Viral Biomass Objective Function (VBOF).

HMPV has two main glycosylated surface proteins:
1. F (Fusion) protein - 3 N-linked glycosylation sites
2. G (Attachment) protein - 5 N-linked + extensive O-linked glycosylation

Key Functions:
- calculate_glycan_stoichiometry: Calculate total glycan requirements for VBOF
- get_f_protein_glycans: Calculate F protein N-linked glycans
- get_g_protein_glycans: Calculate G protein N-linked and O-linked glycans

=============================================================================
SOURCES AND REFERENCES:
=============================================================================

1. F PROTEIN GLYCOSYLATION:
   - Viswanathan et al. (2011) "Effects of N-linked glycosylation of the fusion 
     protein on replication of human metapneumovirus in vitro and in mouse lungs"
     J Gen Virol 92:1580-1584. DOI: 10.1099/vir.0.030049-0
     https://www.researchgate.net/publication/50936404
   
   Key findings:
   - F protein has 3 N-linked glycosylation sites (all utilized)
   - Sites at positions: N57, N172, N353 (or N58, N172, N350 depending on strain)
   - N172 glycan is critical for viral replication
   - Mutations at glycosylation sites impair viral growth

2. G PROTEIN GLYCOSYLATION:
   - Thammawat et al. (2008) "Intracellular Processing, Glycosylation, and Cell 
     Surface Expression of Human Metapneumovirus Attachment Glycoprotein"
     J Virol 82(20):10022-10034. DOI: 10.1128/JVI.01287-06
     https://journals.asm.org/doi/10.1128/jvi.01287-06
     https://pmc.ncbi.nlm.nih.gov/articles/PMC2168831/
   
   Key findings:
   - G protein has 5 potential N-linked glycosylation sites
   - G protein is mucin-like with 30-34% Ser/Thr content
   - Extensive O-linked glycosylation
   - O-glycosylation initiates in trans-Golgi compartment
   - Mature G protein: ~80 kDa vs 27 kDa backbone (glycans add ~53 kDa)

3. GLYCAN STRUCTURE REFERENCES:
   - Varki et al. (2017) "Essentials of Glycobiology" 3rd Edition, NCBI
     https://www.ncbi.nlm.nih.gov/books/NBK310274/
   - High-mannose N-glycan: Man5-9GlcNAc2
   - Complex N-glycan: Core (Man3GlcNAc2) + branches (GlcNAc, Gal, Neu5Ac, Fuc)

=============================================================================
"""

import logging
from typing import Dict, Tuple
from dataclasses import dataclass

from .exceptions import HMPVModelError

logger = logging.getLogger(__name__)


# =============================================================================
# GLYCAN COMPOSITION DATA FROM LITERATURE
# =============================================================================

# F protein glycosylation data
# Sources:
# - Viswanathan et al. (2011) J Gen Virol 92:1580-1584
# - Schowalter et al. (2006) J Virol 80:10931-10941 (PMC1642150)
# - Nature Communications 2017 (cryo-EM prefusion F structure)
#
# KEY FINDINGS FROM MULTIPLE PAPERS:
# - 3 N-linked glycosylation sites: N57, N172, N353 (or N58, N172, N350 depending on strain)
# - ALL THREE SITES ARE UTILIZED (confirmed by Schowalter 2006)
# - Each site affects cleavage and fusion to various degrees
# - N57 and N172 mutations impair viral replication in vitro AND in vivo
# - Dense glycan shield at apex of prefusion F (cryo-EM)
F_PROTEIN_GLYCAN_DATA = {
    'n_linked_sites': 3,  # N57, N172, N353
    'site_positions': [57, 172, 353],  # Asparagine positions (can vary by strain)
    'site_positions_alt': [58, 172, 350],  # Alternative numbering in some strains
    'site_utilization': 1.0,  # All 3 sites are utilized (100%) - confirmed
    'glycan_type': 'complex',  # Complex-type N-glycan on mature F protein
    'critical_sites': [57, 172],  # Critical for replication (Viswanathan 2011)
    'glycan_shield': True,  # Dense glycan shield at apex (cryo-EM)
    
    # Primary source
    'source': 'Viswanathan et al. (2011) J Gen Virol 92:1580-1584',
    'doi': '10.1099/vir.0.030049-0',
    'url': 'https://www.researchgate.net/publication/50936404',
    
    # Additional confirming sources
    'additional_sources': [
        {
            'authors': 'Schowalter et al.',
            'title': 'Characterization of human metapneumovirus F protein-promoted membrane fusion',
            'journal': 'J Virol',
            'year': 2006,
            'volume': '80(22)',
            'pages': '10931-10941',
            'doi': '10.1128/JVI.01287-06',
            'pmc': 'PMC1642150',
            'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC1642150/',
            'key_finding': 'All 3 N-glycosylation sites are utilized'
        },
        {
            'authors': 'Battles et al.',
            'title': 'Structure and immunogenicity of pre-fusion-stabilized human metapneumovirus F glycoprotein',
            'journal': 'Nat Commun',
            'year': 2017,
            'doi': '10.1038/s41467-017-01708-9',
            'url': 'https://www.nature.com/articles/s41467-017-01708-9',
            'key_finding': 'Dense glycan shield at apex of prefusion F'
        }
    ]
}

# G protein glycosylation data
# Sources: 
# - Thammawat et al. (2008) J Virol 82:10022-10034
# - MDPI Viruses review (2014) https://www.mdpi.com/1999-4915/6/8/3019
# - Recent MDPI Disease review (2024)
#
# KEY FINDINGS FROM MULTIPLE PAPERS:
# - N-linked sites: 2-5 (genotype dependent)
#   - Sites at aa 30 and 52 are conserved in genotype A
# - O-linked sites: >60 POTENTIAL sites (higher than previously estimated)
#   - Serine residues: 17-43 (genotype A has more)
#   - Threonine residues: 33-51 (genotype B has more)
#   - Some variants (A2.2.2) have even more due to gene duplications
# - Mature G protein: ~80 kDa vs ~27 kDa backbone
G_PROTEIN_GLYCAN_DATA = {
    'n_linked_sites': 4,  # 2-5 sites (using average of 4)
    'n_linked_sites_range': (2, 5),  # Range across genotypes
    'n_linked_conserved': [30, 52],  # Conserved in genotype A
    'n_linked_utilization': 0.8,  # Estimated ~80% occupancy
    
    # O-linked glycosylation (UPDATED based on literature)
    'o_linked_sites_potential': 60,  # >60 potential sites per MDPI review
    'serine_residues_range': (17, 43),  # Serine count varies by genotype
    'threonine_residues_range': (33, 51),  # Threonine count varies
    'o_linked_sites': 45,  # Conservative estimate: ~45 actual sites
    'o_linked_utilization': 0.65,  # ~65% occupancy (mucin-like)
    
    'ser_thr_percent': 34,  # 30-34% of G protein is Ser/Thr
    'mature_mw_kda': 80,  # Mature glycosylated G protein
    'backbone_mw_kda': 27,  # Unglycosylated backbone
    'glycan_mass_kda': 53,  # ~53 kDa from glycans (mostly O-linked)
    
    # Sources
    'source': 'Thammawat et al. (2008) J Virol 82:10022-10034',
    'doi': '10.1128/JVI.01287-06',
    'url': 'https://journals.asm.org/doi/10.1128/jvi.01287-06',
    'pmc': 'PMC2168831',
    'additional_sources': [
        {
            'title': 'Human metapneumovirus glycoprotein G mediates viral attachment',
            'journal': 'Viruses',
            'year': 2014,
            'url': 'https://www.mdpi.com/1999-4915/6/8/3019'
        }
    ]
}


# =============================================================================
# GLYCAN MONOSACCHARIDE COMPOSITION
# =============================================================================

# Typical N-linked glycan composition (complex type)
# Based on viral glycoprotein glycomics studies
# Source: Varki et al. (2017) Essentials of Glycobiology
N_LINKED_GLYCAN_COMPLEX = {
    'name': 'Complex-type N-glycan (bi-antennary)',
    'composition': {
        'GlcNAc': 4,    # 2 core + 2 antennae
        'Man': 3,       # Core trimannose
        'Gal': 2,       # Terminal galactose on each antenna
        'Neu5Ac': 2,    # Sialic acid caps (common on human glycoproteins)
        'Fuc': 1,       # Core fucose (common)
    },
    'molecular_weight': 2200,  # Average MW in Daltons
    'source': 'Varki et al. (2017) Essentials of Glycobiology, 3rd Ed'
}

# High-mannose N-glycan (often found on viral proteins)
N_LINKED_GLYCAN_HIGH_MANNOSE = {
    'name': 'High-mannose N-glycan (Man9)',
    'composition': {
        'GlcNAc': 2,    # Core chitobiose
        'Man': 9,       # High-mannose
    },
    'molecular_weight': 1900,  # Average MW in Daltons
    'source': 'Varki et al. (2017) Essentials of Glycobiology, 3rd Ed'
}

# O-linked glycan (Core 1, mucin-type)
# Common on mucin-like proteins like HMPV G
O_LINKED_GLYCAN_CORE1 = {
    'name': 'Core 1 O-glycan (T antigen, sialylated)',
    'composition': {
        'GalNAc': 1,    # N-acetylgalactosamine (linked to Ser/Thr)
        'Gal': 1,       # Galactose
        'Neu5Ac': 1,    # Sialic acid cap (common sialyl-T antigen)
    },
    'molecular_weight': 656,  # Average MW in Daltons
    'source': 'Varki et al. (2017) Essentials of Glycobiology, 3rd Ed'
}


# =============================================================================
# BiGG METABOLITE IDs FOR GLYCAN PRECURSORS
# =============================================================================

# Nucleotide sugar donors (activated forms used in glycosylation)
# These are the metabolic precursors consumed during glycan synthesis
GLYCAN_PRECURSOR_BIGG_IDS = {
    # N-linked glycan precursors
    'GlcNAc': 'uacgam_c',     # UDP-N-acetyl-D-glucosamine
    'Man': 'gdpmann_c',       # GDP-mannose
    'Gal': 'udpgal_c',        # UDP-galactose
    'Fuc': 'gdpfuc_c',        # GDP-fucose
    'Neu5Ac': 'cmpacna_c',    # CMP-N-acetylneuraminate (sialic acid)
    
    # O-linked glycan precursors
    'GalNAc': 'udpgalfur_c',  # UDP-N-acetyl-D-galactosamine (approximate)
    
    # Alternative IDs (model-specific)
    'alt_GlcNAc': 'udpacgal_c',
    'alt_Man': 'man_c',
}

# Energy requirements for glycosylation
# Each glycosidic bond requires energy
GLYCOSYLATION_ENERGY = {
    'atp_per_n_glycan': 2,  # ATP for oligosaccharyltransferase
    'atp_per_o_glycan': 1,  # ATP for O-GalNAc transferase
}


@dataclass
class GlycanStoichiometry:
    """Container for glycan stoichiometry data."""
    f_protein_glycans: Dict[str, int]
    g_protein_glycans: Dict[str, int]
    total_glycans: Dict[str, int]
    metadata: Dict


def calculate_f_protein_glycans(
    copy_number: int = 350,
    glycan_type: str = 'complex'
) -> Dict[str, int]:
    """
    Calculate glycan requirements for HMPV F protein.
    
    The F protein has 3 N-linked glycosylation sites that are all utilized.
    
    Parameters:
    -----------
    copy_number : int
        Number of F protein copies per virion (default: 350)
    glycan_type : str
        Type of N-glycan: 'complex' or 'high_mannose' (default: 'complex')
    
    Returns:
    --------
    dict : Monosaccharide requirements (negative = consumption)
    
    Source:
    -------
    Viswanathan et al. (2011) J Gen Virol 92:1580-1584
    - F protein has 3 N-linked glycosylation sites (N57, N172, N353)
    - All 3 sites are utilized
    """
    n_sites = F_PROTEIN_GLYCAN_DATA['n_linked_sites']  # 3 sites
    utilization = F_PROTEIN_GLYCAN_DATA['site_utilization']  # 1.0 (100%)
    
    # Total N-glycans on all F proteins
    total_n_glycans = int(n_sites * utilization * copy_number)
    
    # Select glycan composition
    if glycan_type == 'high_mannose':
        glycan = N_LINKED_GLYCAN_HIGH_MANNOSE
    else:
        glycan = N_LINKED_GLYCAN_COMPLEX
    
    # Calculate monosaccharide requirements
    glycan_stoichiometry = {}
    for monosaccharide, count_per_glycan in glycan['composition'].items():
        total_count = count_per_glycan * total_n_glycans
        glycan_stoichiometry[monosaccharide] = -total_count  # Negative = consumption
    
    logger.info(f"F protein glycan stoichiometry calculated:")
    logger.info(f"  Copy number: {copy_number}")
    logger.info(f"  N-linked sites: {n_sites}")
    logger.info(f"  Total N-glycans: {total_n_glycans}")
    logger.info(f"  Glycan type: {glycan_type}")
    
    return glycan_stoichiometry


def calculate_g_protein_glycans(
    copy_number: int = 250
) -> Dict[str, int]:
    """
    Calculate glycan requirements for HMPV G protein.
    
    The G protein is heavily glycosylated with both N-linked and O-linked glycans.
    It has a mucin-like structure with 30-34% Ser/Thr content.
    
    Parameters:
    -----------
    copy_number : int
        Number of G protein copies per virion (default: 250)
    
    Returns:
    --------
    dict : Monosaccharide requirements (negative = consumption)
    
    Sources:
    --------
    1. Thammawat et al. (2008) J Virol 82:10022-10034 (PMC2168831)
       - G protein has 5 potential N-linked glycosylation sites
       - Extensive O-linked glycosylation (mucin-like, 30-34% Ser/Thr)
       - Mature G protein is ~80 kDa vs ~27 kDa backbone
    
    2. MDPI Viruses (2014) https://www.mdpi.com/1999-4915/6/8/3019
       - N-linked sites: 2-5 (genotype dependent)
       - O-linked sites: >60 POTENTIAL sites
       - Serine residues: 17-43
       - Threonine residues: 33-51
    
    3. Schowalter et al. (2006) J Virol (PMC1642150)
       - Additional confirmation of G protein glycosylation patterns
    """
    # N-linked glycans (2-5 sites depending on genotype, using average of 4)
    n_sites = G_PROTEIN_GLYCAN_DATA['n_linked_sites']  # 4 sites (average)
    n_utilization = G_PROTEIN_GLYCAN_DATA['n_linked_utilization']  # ~80%
    total_n_glycans = int(n_sites * n_utilization * copy_number)
    
    # O-linked glycans (mucin-like region) - UPDATED based on MDPI review
    # Literature indicates >60 potential sites, we use 45 as conservative estimate
    o_sites = G_PROTEIN_GLYCAN_DATA['o_linked_sites']  # ~45 sites (updated from 25)
    o_utilization = G_PROTEIN_GLYCAN_DATA['o_linked_utilization']  # ~65%
    total_o_glycans = int(o_sites * o_utilization * copy_number)
    
    # Calculate monosaccharide requirements
    glycan_stoichiometry = {}
    
    # N-linked glycans (using complex type)
    for monosaccharide, count_per_glycan in N_LINKED_GLYCAN_COMPLEX['composition'].items():
        total_count = count_per_glycan * total_n_glycans
        glycan_stoichiometry[monosaccharide] = glycan_stoichiometry.get(monosaccharide, 0) - total_count
    
    # O-linked glycans (using Core 1 type)
    for monosaccharide, count_per_glycan in O_LINKED_GLYCAN_CORE1['composition'].items():
        total_count = count_per_glycan * total_o_glycans
        glycan_stoichiometry[monosaccharide] = glycan_stoichiometry.get(monosaccharide, 0) - total_count
    
    logger.info(f"G protein glycan stoichiometry calculated:")
    logger.info(f"  Copy number: {copy_number}")
    logger.info(f"  N-linked sites: {n_sites}, Total N-glycans: {total_n_glycans}")
    logger.info(f"  O-linked sites: ~{o_sites}, Total O-glycans: {total_o_glycans}")
    
    return glycan_stoichiometry


def calculate_glycan_stoichiometry(
    f_copy_number: int = 350,
    g_copy_number: int = 250,
    include_sh_glycans: bool = False
) -> Tuple[Dict[str, float], Dict]:
    """
    Calculate total glycan requirements for HMPV VBOF.
    
    Combines glycan requirements from all glycosylated proteins:
    - F protein: 3 N-linked glycosylation sites
    - G protein: 5 N-linked + ~25 O-linked glycosylation sites
    - SH protein: Minimal glycosylation (optional)
    
    Parameters:
    -----------
    f_copy_number : int
        Number of F protein copies per virion
    g_copy_number : int
        Number of G protein copies per virion
    include_sh_glycans : bool
        Whether to include SH protein glycosylation (default: False)
    
    Returns:
    --------
    tuple : (stoichiometry_dict, metadata_dict)
        - stoichiometry_dict: BiGG metabolite IDs with coefficients
        - metadata_dict: Summary information
    
    Sources:
    --------
    1. Viswanathan et al. (2011) J Gen Virol 92:1580-1584
       DOI: 10.1099/vir.0.030049-0
       https://www.researchgate.net/publication/50936404
    
    2. Thammawat et al. (2008) J Virol 82:10022-10034
       DOI: 10.1128/JVI.01287-06
       https://journals.asm.org/doi/10.1128/jvi.01287-06
    """
    logger.info("Calculating glycan stoichiometry for HMPV VBOF...")
    
    # Calculate glycans for each protein
    f_glycans = calculate_f_protein_glycans(f_copy_number)
    g_glycans = calculate_g_protein_glycans(g_copy_number)
    
    # Combine all glycan requirements
    total_monosaccharides = {}
    for monosaccharide, count in f_glycans.items():
        total_monosaccharides[monosaccharide] = total_monosaccharides.get(monosaccharide, 0) + count
    for monosaccharide, count in g_glycans.items():
        total_monosaccharides[monosaccharide] = total_monosaccharides.get(monosaccharide, 0) + count
    
    # Convert to BiGG IDs (nucleotide sugar donors)
    bigg_stoichiometry = {}
    for monosaccharide, count in total_monosaccharides.items():
        bigg_id = GLYCAN_PRECURSOR_BIGG_IDS.get(monosaccharide)
        if bigg_id:
            bigg_stoichiometry[bigg_id] = count
        else:
            logger.warning(f"No BiGG ID mapping for monosaccharide: {monosaccharide}")
            # Use placeholder
            bigg_stoichiometry[f'{monosaccharide.lower()}_c'] = count
    
    # Add energy requirements for glycosylation
    total_n_glycans = (F_PROTEIN_GLYCAN_DATA['n_linked_sites'] * f_copy_number +
                       int(G_PROTEIN_GLYCAN_DATA['n_linked_sites'] * 
                           G_PROTEIN_GLYCAN_DATA['n_linked_utilization'] * g_copy_number))
    total_o_glycans = int(G_PROTEIN_GLYCAN_DATA['o_linked_sites'] * 
                         G_PROTEIN_GLYCAN_DATA['o_linked_utilization'] * g_copy_number)
    
    glycosylation_atp = (total_n_glycans * GLYCOSYLATION_ENERGY['atp_per_n_glycan'] +
                         total_o_glycans * GLYCOSYLATION_ENERGY['atp_per_o_glycan'])
    
    # Create metadata
    metadata = {
        'f_protein': {
            'copy_number': f_copy_number,
            'n_linked_sites': F_PROTEIN_GLYCAN_DATA['n_linked_sites'],
            'total_n_glycans': F_PROTEIN_GLYCAN_DATA['n_linked_sites'] * f_copy_number,
            'source': F_PROTEIN_GLYCAN_DATA['source'],
            'doi': F_PROTEIN_GLYCAN_DATA['doi']
        },
        'g_protein': {
            'copy_number': g_copy_number,
            'n_linked_sites': G_PROTEIN_GLYCAN_DATA['n_linked_sites'],
            'o_linked_sites_estimated': G_PROTEIN_GLYCAN_DATA['o_linked_sites'],
            'total_n_glycans': total_n_glycans - F_PROTEIN_GLYCAN_DATA['n_linked_sites'] * f_copy_number,
            'total_o_glycans': total_o_glycans,
            'source': G_PROTEIN_GLYCAN_DATA['source'],
            'doi': G_PROTEIN_GLYCAN_DATA['doi']
        },
        'total_monosaccharides': {k: abs(v) for k, v in total_monosaccharides.items()},
        'glycosylation_atp': glycosylation_atp,
        'total_glycan_sites': total_n_glycans + total_o_glycans
    }
    
    logger.info(f"Glycan stoichiometry calculation complete:")
    logger.info(f"  Total N-glycan sites: {total_n_glycans}")
    logger.info(f"  Total O-glycan sites: {total_o_glycans}")
    logger.info(f"  Glycosylation ATP: {glycosylation_atp}")
    for mono, count in total_monosaccharides.items():
        logger.info(f"  {mono}: {abs(count)}")
    
    return bigg_stoichiometry, metadata


def get_glycan_summary() -> str:
    """
    Get a summary of glycan data sources and calculations.
    
    Returns:
    --------
    str : Formatted summary string
    """
    summary = """
=============================================================================
HMPV GLYCAN DATA SUMMARY
=============================================================================

F PROTEIN GLYCOSYLATION:
------------------------
Primary Source: Viswanathan et al. (2011) J Gen Virol 92:1580-1584
DOI: 10.1099/vir.0.030049-0
URL: https://www.researchgate.net/publication/50936404

Confirming Source: Schowalter et al. (2006) J Virol 80:10931-10941
DOI: 10.1128/JVI.01287-06
PMC: PMC1642150
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC1642150/

Structural Source: Battles et al. (2017) Nat Commun 8:1528
DOI: 10.1038/s41467-017-01708-9
URL: https://www.nature.com/articles/s41467-017-01708-9

Key findings:
- 3 N-linked glycosylation sites: N57, N172, N353
- ALL 3 SITES ARE UTILIZED (100% occupancy) - Schowalter 2006
- N57 and N172 glycans critical for viral replication in vitro AND in vivo
- Dense glycan shield at apex of prefusion F (cryo-EM - Battles 2017)
- Each glycan affects cleavage and fusion to various degrees

G PROTEIN GLYCOSYLATION:
------------------------
Primary Source: Thammawat et al. (2008) J Virol 82:10022-10034
DOI: 10.1128/JVI.01287-06
PMC: PMC2168831
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC2168831/

Additional Source: MDPI Viruses Review (2014)
URL: https://www.mdpi.com/1999-4915/6/8/3019

Key findings:
- N-linked glycosylation sites: 2-5 (genotype dependent)
  - Conserved sites at aa 30 and 52 in genotype A
- O-linked glycosylation: >60 POTENTIAL SITES
  - Serine residues: 17-43 (genotype A has more)
  - Threonine residues: 33-51 (genotype B has more)
- Mucin-like structure with 30-34% Ser/Thr content
- Mature G protein: ~80 kDa (backbone ~27 kDa + ~53 kDa glycans)
- O-glycosylation initiates in trans-Golgi compartment

GLYCAN COMPOSITION:
-------------------
N-linked (complex type): GlcNAc4-Man3-Gal2-Neu5Ac2-Fuc1 (~2200 Da)
O-linked (Core 1 sialyl-T): GalNAc1-Gal1-Neu5Ac1 (~656 Da)

=============================================================================
"""
    return summary

