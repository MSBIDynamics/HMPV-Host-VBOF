"""
Glycan Analyzer Module
======================

This module calculates glycan (carbohydrate) requirements for HMPV glycoproteins
for inclusion in the Viral Biomass Objective Function (VBOF).


Key Functions:
- calculate_glycan_stoichiometry: Calculate total glycan requirements for VBOF
- get_f_protein_glycans: Calculate F protein N-linked glycans
- get_g_protein_glycans: Calculate G protein N-linked and O-linked glycans

Author: Syed Mushahid Hussain

=============================================================================
SOURCES AND REFERENCES:
=============================================================================
    HMPV has two main glycosylated surface proteins:
    1. F (Fusion) protein 
    2. G (Attachment) protein

  
   - F protein has 3 N-linked glycosylation sites (all utilized)
   - Sites at positions: N57, N172, N353 (or N58, N172, N350 depending on strain)
   Sources:  https://journals.asm.org/doi/10.1128/jvi.01287-06
   

   - sG protein has 3 to 6 potential N-linked glycosylation sites and 26 O-linked sites, but utilization is variable and likely incomplete.
    Sources:  https://journals.asm.org/doi/10.1128/jvi.01469-07
    
   - G protein has 5 potential N-linked sites and 59 O-linked sites.
   Sources:  https://journals.asm.org/doi/10.1128/jvi.01726-14
  
   - G protein is mucin-like with 30-34% Ser/Thr content
   - Extensive O-linked glycosylation
   Sources:  https://journals.asm.org/doi/10.1128/jvi.01469-07

   GLYCAN STRUCTURE REFERENCES:
   - Varki et al. (2017) "Essentials of Glycobiology" 3rd Edition, NCBI
     https://www.ncbi.nlm.nih.gov/books/NBK310274/
   - High-mannose N-glycan: Man5-9GlcNAc2
   - Complex N-glycan: Core (Man3GlcNAc2) + branches (GlcNAc, Gal, Neu5Ac, Fuc)

=============================================================================
"""

import logging
from typing import Dict, Tuple
from dataclasses import dataclass



from .config import (
    F_PROTEIN_N_LINKED_SITES,
    F_PROTEIN_SITE_UTILIZATION,
    G_PROTEIN_N_LINKED_SITES,
    G_PROTEIN_N_LINKED_UTILIZATION,
    G_PROTEIN_O_LINKED_SITES,
    G_PROTEIN_O_LINKED_UTILIZATION,
    N_GLYCAN_COMPLEX_COMPOSITION,
    N_LINKED_GLYCAN_HIGH_MANNOSE,
    O_GLYCAN_CORE1_COMPOSITION,
    GLYCOSYLATION_ATP_PER_N_GLYCAN,
    GLYCOSYLATION_ATP_PER_O_GLYCAN,
    GLYCAN_PRECURSOR_BIGG_IDS,
)

logger = logging.getLogger(__name__)


@dataclass
class GlycanStoichiometry:
    """Container for glycan stoichiometry data."""
    f_protein_glycans: Dict[str, int]
    g_protein_glycans: Dict[str, int]
    total_glycans: Dict[str, int]
    metadata: Dict


def calculate_f_protein_glycans(
    copy_number: int = 693,
    glycan_type: str = 'complex'
) -> Dict[str, int]:
    """
    Calculate glycan requirements for HMPV F protein.
    
    The F protein has 3 N-linked glycosylation sites that are all utilized.
    
    Parameters:
    -----------
    copy_number : int
        Number of F protein copies per virion (default: 693)
    glycan_type : str
        Type of N-glycan: 'complex' or 'high_mannose' (default: 'complex')
    
    Returns:
    --------
    dict : Monosaccharide requirements (negative = consumption)
    

    """
    n_sites = F_PROTEIN_N_LINKED_SITES        
    utilization = F_PROTEIN_SITE_UTILIZATION  

    # Total N-glycans on all F proteins
    total_n_glycans = int(n_sites * utilization * copy_number)

    # Select glycan composition
    if glycan_type == 'high_mannose':
        glycan_composition = N_LINKED_GLYCAN_HIGH_MANNOSE
    else:
        glycan_composition = N_GLYCAN_COMPLEX_COMPOSITION  

    # Calculate monosaccharide requirements
    glycan_stoichiometry = {}
    for monosaccharide, count_per_glycan in glycan_composition.items():
        total_count = count_per_glycan * total_n_glycans
        glycan_stoichiometry[monosaccharide] = -total_count  # Negative = consumption
    
    logger.info(f"F protein glycan stoichiometry calculated:")
    logger.info(f"  Copy number: {copy_number}")
    logger.info(f"  N-linked sites: {n_sites}")
    logger.info(f"  Total N-glycans: {total_n_glycans}")
    logger.info(f"  Glycan type: {glycan_type}")
    
    return glycan_stoichiometry


def calculate_g_protein_glycans(
    copy_number: int = 139
) -> Dict[str, int]:
    """
    Calculate glycan requirements for HMPV G protein.
    
    The G protein is heavily glycosylated with both N-linked and O-linked glycans.
    It has a mucin-like structure with 30-34% Ser/Thr content.
    
    Parameters:
    -----------
    copy_number : int
        Number of G protein copies per virion (default: 139)
    
    Returns:
    --------
    dict : Monosaccharide requirements (negative = consumption)
    
    
    """
    # N-linked glycans — controlled via config
    n_sites = G_PROTEIN_N_LINKED_SITES           # from config (default 4)
    n_utilization = G_PROTEIN_N_LINKED_UTILIZATION  # from config (default 0.8)
    total_n_glycans = int(n_sites * n_utilization * copy_number)

    # O-linked glycans — controlled via config
    o_sites = G_PROTEIN_O_LINKED_SITES           # from config (default 26)
    o_utilization = G_PROTEIN_O_LINKED_UTILIZATION  # from config (default 0.50)
    total_o_glycans = int(o_sites * o_utilization * copy_number)

    # Calculate monosaccharide requirements
    glycan_stoichiometry = {}

    # N-linked glycans (complex type) — composition from config
    for monosaccharide, count_per_glycan in N_GLYCAN_COMPLEX_COMPOSITION.items():
        total_count = count_per_glycan * total_n_glycans
        glycan_stoichiometry[monosaccharide] = glycan_stoichiometry.get(monosaccharide, 0) - total_count

    # O-linked glycans (Core 1 type) — composition from config
    for monosaccharide, count_per_glycan in O_GLYCAN_CORE1_COMPOSITION.items():
        total_count = count_per_glycan * total_o_glycans
        glycan_stoichiometry[monosaccharide] = glycan_stoichiometry.get(monosaccharide, 0) - total_count
    
    logger.info(f"G protein glycan stoichiometry calculated:")
    logger.info(f"  Copy number: {copy_number}")
    logger.info(f"  N-linked sites: {n_sites}, Total N-glycans: {total_n_glycans}")
    logger.info(f"  O-linked sites: ~{o_sites}, Total O-glycans: {total_o_glycans}")
    
    return glycan_stoichiometry


def calculate_glycosylation_atp(
    f_copy_number: int = 693,
    g_copy_number: int = 139,
) -> float:
    """
    ATP cost for N- and O-linked glycosylation of F and G (transfer machinery).

    Uses the same site counts as calculate_glycan_stoichiometry; rates from config.
    """
    total_n_glycans = (
        F_PROTEIN_N_LINKED_SITES * f_copy_number
        + int(G_PROTEIN_N_LINKED_SITES * G_PROTEIN_N_LINKED_UTILIZATION * g_copy_number)
    )
    total_o_glycans = int(
        G_PROTEIN_O_LINKED_SITES * G_PROTEIN_O_LINKED_UTILIZATION * g_copy_number
    )
    return float(
        total_n_glycans * GLYCOSYLATION_ATP_PER_N_GLYCAN
        + total_o_glycans * GLYCOSYLATION_ATP_PER_O_GLYCAN
    )


def calculate_glycan_stoichiometry(
    f_copy_number: int = 693,
    g_copy_number: int = 139,
    include_sh_glycans: bool = False
) -> Tuple[Dict[str, float], Dict]:
    """
    Calculate total glycan requirements for HMPV VBOF.
    This function combines F and G protein glycan requirements.
    
    Parameters:
    -----------
    f_copy_number : int
        Number of F protein copies per virion
    g_copy_number : int
        Number of G protein copies per virion
    include_sh_glycans : bool
        Whether to include SH protein glycosylation (default: False) Very low number and likely minimal contribution, so excluded by default.
    
    Returns:
    --------
    tuple : (stoichiometry_dict, metadata_dict)
        - stoichiometry_dict: BiGG metabolite IDs with coefficients
        - metadata_dict: Summary information
    

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
    
    # Energy for glycosylation (also applied in vbof_builder.calculate_energy_requirements)
    total_n_glycans = (
        F_PROTEIN_N_LINKED_SITES * f_copy_number
        + int(G_PROTEIN_N_LINKED_SITES * G_PROTEIN_N_LINKED_UTILIZATION * g_copy_number)
    )
    total_o_glycans = int(
        G_PROTEIN_O_LINKED_SITES * G_PROTEIN_O_LINKED_UTILIZATION * g_copy_number
    )
    glycosylation_atp = calculate_glycosylation_atp(f_copy_number, g_copy_number)
    
    # Create metadata
    metadata = {
        'f_protein': {
            'copy_number': f_copy_number,
            'n_linked_sites': F_PROTEIN_N_LINKED_SITES,
            'total_n_glycans': F_PROTEIN_N_LINKED_SITES * f_copy_number,
            
        },
        'g_protein': {
            'copy_number': g_copy_number,
            'n_linked_sites': G_PROTEIN_N_LINKED_SITES,
            'o_linked_sites_estimated': G_PROTEIN_O_LINKED_SITES,
            'total_n_glycans': total_n_glycans - F_PROTEIN_N_LINKED_SITES * f_copy_number,
            'total_o_glycans': total_o_glycans,
           
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

