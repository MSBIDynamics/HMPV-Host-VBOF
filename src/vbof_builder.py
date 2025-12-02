"""
VBOF Builder Module
===================

This module assembles all components into a complete Viral Biomass Objective Function (VBOF)
for HMPV and integrates it with the host metabolic model.

Key Functions:
- calculate_energy_requirements: Calculate ATP/GTP costs
- calculate_lipid_stoichiometry: Calculate envelope lipid requirements
- build_vbof: Assemble complete VBOF stoichiometry
- create_vbof_reaction: Create COBRApy reaction object

Author: Syed Mushahid Hussain
"""

import logging
import math
from typing import Dict, Optional, Tuple
from dataclasses import dataclass

from .exceptions import VBOFConstructionError, EnergyDataMissingError
from .glycan_analyzer import calculate_glycan_stoichiometry

logger = logging.getLogger(__name__)


@dataclass
class VBOFComponents:
    """Container for all VBOF components."""
    genome_stoichiometry: Dict[str, float]
    protein_stoichiometry: Dict[str, float]
    energy_stoichiometry: Dict[str, float]
    lipid_stoichiometry: Dict[str, float]
    glycan_stoichiometry: Dict[str, float]  # NEW: Glycan requirements
    combined_stoichiometry: Dict[str, float]
    metadata: Dict


# =============================================================================
# ENERGY COSTS - Based on biochemistry literature
# =============================================================================
#
# SOURCES AND REFERENCES:
#
# 1. RNA POLYMERIZATION ENERGY:
#    - Berg, Tymoczko, Stryer (2015) "Biochemistry" 8th Edition
#      Chapter 29: RNA Synthesis and Processing
#      https://www.ncbi.nlm.nih.gov/books/NBK22445/
#    - NTP → NMP + PPi (equivalent to 2 ATP per nucleotide)
#
# 2. PROTEIN SYNTHESIS ENERGY:
#    - Alberts et al. (2015) "Molecular Biology of the Cell" 6th Edition
#      Chapter 6: How Cells Read the Genome
#      https://www.ncbi.nlm.nih.gov/books/NBK26887/
#    - tRNA charging: AA + ATP → AA-AMP → AA-tRNA (2 ATP equiv)
#    - Elongation: EF-Tu·GTP + EF-G·GTP (2 GTP per amino acid)
#
# 3. VBOF METHODOLOGY VALIDATION:
#    - Aller et al. (2018) "Integrated human-virus metabolic stoichiometric modelling"
#      https://pmc.ncbi.nlm.nih.gov/articles/PMC6170780/
#
# 4. SARS-CoV-2 VBOF ENERGY CALCULATIONS:
#    - "Integration of metabolic networks and gene expression"
#      https://academic.oup.com/bioinformatics/article/21/18/3645/202235?login=false
#
ENERGY_COSTS = {
    'rna_polymerization_atp_per_nt': 2,  # 2 ATP equivalents per nucleotide 
    'trna_charging_atp_per_aa': 2,        # Amino acid activation (Alberts 2015)
    'elongation_gtp_per_aa': 2,           # EF-Tu and EF-G (Alberts 2015)
    'initiation_gtp': 1,                  # Initiation factor IF-2
    'termination_gtp': 1,                 # Release factor RF-3
}


def calculate_energy_requirements(
    genome_length: int,
    total_amino_acids: int,
    num_proteins: int = 9
) -> Dict[str, float]:
    """
    Calculate ATP and GTP requirements for viral replication.
    
    Energy costs include:
    1. RNA synthesis: NTP → NMP + PPi (2 ATP equiv per nt)
    2. tRNA charging: AA + ATP → AA-tRNA + AMP + PPi (2 ATP per aa)
    3. Translation elongation: 2 GTP per amino acid
    4. Initiation/termination: additional GTP costs
    
    Parameters:
    -----------
    genome_length : int
        Length of viral genome in nucleotides
    total_amino_acids : int
        Total amino acids in all proteins × copy numbers
    num_proteins : int
        Number of distinct proteins (for init/term calculation)
    
    Returns:
    --------
    dict : Energy stoichiometry (negative = consumption)
    
    Note:
    -----
    These are minimum estimates. Actual costs may be higher due to
    proofreading, chaperones, and other processes.
    
    Source:
    -------
    Standard biochemistry values, validated in:
    - Aller et al. (2018) J R Soc Interface
    - Thiele et al. (2020) Bioinformatics
    """
    if genome_length <= 0:
        raise EnergyDataMissingError("Invalid genome length for energy calculation")
    if total_amino_acids <= 0:
        raise EnergyDataMissingError("Invalid amino acid count for energy calculation")
    
    # RNA synthesis energy
    # Note: This is for the packaged genome, not replication intermediates
    rna_energy_atp = genome_length * ENERGY_COSTS['rna_polymerization_atp_per_nt']
    
    # Protein synthesis energy
    # tRNA charging: 2 ATP per amino acid (AA + ATP → AA-AMP → AA-tRNA)
    protein_atp = total_amino_acids * ENERGY_COSTS['trna_charging_atp_per_aa']
    
    # Elongation: 2 GTP per amino acid (EF-Tu + EF-G)
    protein_gtp = total_amino_acids * ENERGY_COSTS['elongation_gtp_per_aa']
    
    # Initiation and termination (approximate, per protein molecule)
    # This is a simplification - actual value depends on number of protein copies
    init_term_gtp = num_proteins * (ENERGY_COSTS['initiation_gtp'] + ENERGY_COSTS['termination_gtp'])
    
    # Total energy
    total_atp = rna_energy_atp + protein_atp
    total_gtp = protein_gtp + init_term_gtp
    
    # Stoichiometry (negative = consumption, positive = production)
    stoichiometry = {
        # ATP consumption and hydrolysis products
        'atp_c': -total_atp,
        'adp_c': total_atp,      # ADP produced
        'amp_c': protein_atp,    # AMP from tRNA charging
        
        # GTP consumption and hydrolysis products
        # Note: GTP for genome is already in genome_stoichiometry
        # This is additional GTP for translation
        'gtp_c': -total_gtp,     # Additional GTP for translation
        'gdp_c': total_gtp,      # GDP produced
        
        # Phosphate
        'pi_c': total_atp + total_gtp,  # Inorganic phosphate
        
        # Water consumption for hydrolysis
        'h2o_c': -(total_atp + total_gtp),
        
        # Protons (simplified - actual stoichiometry is complex)
        'h_c': total_atp + total_gtp,
    }
    
    logger.info(f"Energy requirements calculated:")
    logger.info(f"  RNA synthesis ATP: {rna_energy_atp}")
    logger.info(f"  Protein synthesis ATP: {protein_atp}")
    logger.info(f"  Translation GTP: {protein_gtp}")
    logger.info(f"  Total ATP: {total_atp}")
    logger.info(f"  Total GTP: {total_gtp}")
    
    return stoichiometry


def calculate_lipid_stoichiometry(
    virion_diameter_nm: float = 150.0,
    include_lipids: bool = True
) -> Dict[str, float]:
    """
    Calculate lipid requirements for the viral envelope.
    
    For HMPV (enveloped paramyxovirus):
    - Envelope derived from host plasma membrane
    - Typical composition similar to other paramyxoviruses
    
    Parameters:
    -----------
    virion_diameter_nm : float
        Average virion diameter in nanometers (default: 150 nm for HMPV)
    include_lipids : bool
        Whether to include lipid requirements (default: True)
    
    Returns:
    --------
    dict : Lipid stoichiometry (negative = consumption)
    
    Note:
    -----
    HMPV can be pleomorphic (variable shape), so this is an approximation
    based on spherical geometry.
    
    Lipid composition based on paramyxovirus literature:
    - Phosphatidylcholine (PC): ~45%
    - Phosphatidylethanolamine (PE): ~25%
    - Phosphatidylserine (PS): ~5%
    - Sphingomyelin (SM): ~15%
    - Cholesterol: ~10%
    
    Sources:
    --------
    1. Jinjun Shan. (2028) "High-resolution lipidomics reveals dysregulation of lipid metabolism"
       https://pubs.rsc.org/en/content/articlelanding/2018/ra/c8ra05640d
    
    2. Brügger et al. (2006) "The HIV lipidome"
       https://pmc.ncbi.nlm.nih.gov/articles/PMC1413831/
    
    3. Nagle & Tristram-Nagle (2000) "Structure of lipid bilayers"
       Biochim Biophys Acta 1469(3):159-195. DOI: 10.1016/S0304-4157(00)00016-2
       https://www.sciencedirect.com/science/article/pii/S0304415700000162
       - Lipid packing density: ~0.65 nm² per lipid molecule
    
    4. Aller et al. (2018) "Integrated human-virus metabolic stoichiometric modelling"
       https://pmc.ncbi.nlm.nih.gov/articles/PMC6170780/
       - VBOF methodology for enveloped viruses
    """
    if not include_lipids:
        logger.info("Lipid stoichiometry skipped (include_lipids=False)")
        return {}
    
    # Calculate surface area (spherical approximation)
    radius_nm = virion_diameter_nm / 2
    surface_area_nm2 = 4 * math.pi * radius_nm ** 2
    
    # Lipid packing density
    # Typical phospholipid area: 0.6-0.7 nm² per lipid
    lipid_area_nm2 = 0.65
    
    # Total lipids in bilayer (2 leaflets)
    total_lipids = int((surface_area_nm2 / lipid_area_nm2) * 2)
    
    # Lipid composition fractions (based on paramyxovirus data)
    # Note: Using generic lipid BiGG IDs - may need model-specific mapping
    lipid_fractions = {
        'pc_hs_c': 0.45,        # Phosphatidylcholine
        'pe_hs_c': 0.25,        # Phosphatidylethanolamine  
        'ps_hs_c': 0.05,        # Phosphatidylserine
        'sphmyln_hs_c': 0.15,   # Sphingomyelin
        'chsterol_c': 0.10,     # Cholesterol
    }
    
    # Calculate stoichiometry
    stoichiometry = {}
    for lipid_id, fraction in lipid_fractions.items():
        count = int(total_lipids * fraction)
        if count > 0:
            stoichiometry[lipid_id] = -count
    
    logger.info(f"Lipid stoichiometry calculated:")
    logger.info(f"  Virion diameter: {virion_diameter_nm} nm")
    logger.info(f"  Surface area: {surface_area_nm2:.0f} nm²")
    logger.info(f"  Total lipids: {total_lipids}")
    for lipid_id, count in stoichiometry.items():
        logger.info(f"  {lipid_id}: {count}")
    
    return stoichiometry


def build_vbof(
    genome_stoichiometry: Dict[str, float],
    protein_stoichiometry: Dict[str, float],
    genome_length: int,
    total_amino_acids: int,
    num_proteins: int = 9,
    virion_diameter_nm: float = 150.0,
    include_lipids: bool = True,
    include_energy: bool = True,
    include_glycans: bool = True,
    f_copy_number: int = 350,
    g_copy_number: int = 250
) -> VBOFComponents:
    """
    Build complete VBOF by combining all components.
    
    Parameters:
    -----------
    genome_stoichiometry : dict
        Nucleotide requirements from genome_analyzer
    protein_stoichiometry : dict
        Amino acid requirements from protein_analyzer
    genome_length : int
        Genome length for energy calculations
    total_amino_acids : int
        Total amino acids for energy calculations
    num_proteins : int
        Number of proteins (default: 9 for HMPV)
    virion_diameter_nm : float
        Virion diameter for lipid calculations
    include_lipids : bool
        Whether to include lipid envelope
    include_energy : bool
        Whether to include energy requirements
    include_glycans : bool
        Whether to include glycan requirements (default: True)
    f_copy_number : int
        F protein copy number for glycan calculation (default: 350)
    g_copy_number : int
        G protein copy number for glycan calculation (default: 250)
    
    Returns:
    --------
    VBOFComponents : Complete VBOF with all components
    
    Raises:
    -------
    VBOFConstructionError : If VBOF cannot be assembled
    
    Glycan Sources:
    ---------------
    1. F protein: Viswanathan et al. (2011) J Gen Virol 92:1580-1584
       - 3 N-linked glycosylation sites (N57, N172, N353)
    2. G protein: Thammawat et al. (2008) J Virol 82:10022-10034
       - 5 N-linked + extensive O-linked glycosylation
    """
    logger.info("Building VBOF...")
    
    # Validate inputs
    if not genome_stoichiometry:
        raise VBOFConstructionError("Genome stoichiometry is empty")
    if not protein_stoichiometry:
        raise VBOFConstructionError("Protein stoichiometry is empty")
    
    # Calculate energy requirements
    if include_energy:
        energy_stoichiometry = calculate_energy_requirements(
            genome_length=genome_length,
            total_amino_acids=total_amino_acids,
            num_proteins=num_proteins
        )
    else:
        energy_stoichiometry = {}
        logger.info("Energy stoichiometry skipped (include_energy=False)")
    
    # Calculate lipid requirements
    lipid_stoichiometry = calculate_lipid_stoichiometry(
        virion_diameter_nm=virion_diameter_nm,
        include_lipids=include_lipids
    )
    
    # Calculate glycan requirements (NEW)
    if include_glycans:
        glycan_stoichiometry, glycan_metadata = calculate_glycan_stoichiometry(
            f_copy_number=f_copy_number,
            g_copy_number=g_copy_number
        )
        logger.info(f"Glycan stoichiometry calculated:")
        logger.info(f"  F protein glycans: {glycan_metadata['f_protein']['total_n_glycans']} N-glycans")
        logger.info(f"  G protein glycans: {glycan_metadata['g_protein']['total_n_glycans']} N-glycans, "
                   f"{glycan_metadata['g_protein']['total_o_glycans']} O-glycans")
    else:
        glycan_stoichiometry = {}
        glycan_metadata = {}
        logger.info("Glycan stoichiometry skipped (include_glycans=False)")
    
    # Combine all stoichiometries
    combined = {}
    
    # Add genome stoichiometry
    for met_id, coef in genome_stoichiometry.items():
        combined[met_id] = combined.get(met_id, 0) + coef
    
    # Add protein stoichiometry
    for met_id, coef in protein_stoichiometry.items():
        combined[met_id] = combined.get(met_id, 0) + coef
    
    # Add energy stoichiometry
    for met_id, coef in energy_stoichiometry.items():
        combined[met_id] = combined.get(met_id, 0) + coef
    
    # Add lipid stoichiometry
    for met_id, coef in lipid_stoichiometry.items():
        combined[met_id] = combined.get(met_id, 0) + coef
    
    # Add glycan stoichiometry (NEW)
    for met_id, coef in glycan_stoichiometry.items():
        combined[met_id] = combined.get(met_id, 0) + coef
    
    # Remove zero coefficients
    combined = {k: v for k, v in combined.items() if v != 0}
    
    # Create metadata
    metadata = {
        'genome_length': genome_length,
        'total_amino_acids': total_amino_acids,
        'num_proteins': num_proteins,
        'virion_diameter_nm': virion_diameter_nm,
        'include_lipids': include_lipids,
        'include_energy': include_energy,
        'include_glycans': include_glycans,
        'total_metabolites': len(combined),
        'total_consumed': sum(1 for v in combined.values() if v < 0),
        'total_produced': sum(1 for v in combined.values() if v > 0),
    }
    
    # Add glycan metadata if available
    if include_glycans and glycan_metadata:
        metadata['glycan_data'] = glycan_metadata
    
    vbof = VBOFComponents(
        genome_stoichiometry=genome_stoichiometry,
        protein_stoichiometry=protein_stoichiometry,
        energy_stoichiometry=energy_stoichiometry,
        lipid_stoichiometry=lipid_stoichiometry,
        glycan_stoichiometry=glycan_stoichiometry,
        combined_stoichiometry=combined,
        metadata=metadata
    )
    
    logger.info(f"VBOF built successfully:")
    logger.info(f"  Total metabolites: {metadata['total_metabolites']}")
    logger.info(f"  Consumed: {metadata['total_consumed']}")
    logger.info(f"  Produced: {metadata['total_produced']}")
    
    return vbof


def get_vbof_summary(vbof: VBOFComponents) -> Dict:
    """
    Get a summary of the VBOF.
    
    Parameters:
    -----------
    vbof : VBOFComponents
        Complete VBOF
    
    Returns:
    --------
    dict : Summary statistics
    """
    consumed = {k: v for k, v in vbof.combined_stoichiometry.items() if v < 0}
    produced = {k: v for k, v in vbof.combined_stoichiometry.items() if v > 0}
    
    # Categorize metabolites
    nucleotides = {k: v for k, v in consumed.items() if k.endswith('tp_c')}
    amino_acids = {k: v for k, v in consumed.items() if '__L_c' in k or k == 'gly_c'}
    lipids = {k: v for k, v in consumed.items() if 'hs_c' in k or 'chsterol' in k}
    energy = {k: v for k, v in consumed.items() if k in ['atp_c', 'gtp_c'] and k not in nucleotides}
    
    summary = {
        'metadata': vbof.metadata,
        'categories': {
            'nucleotides': {
                'count': len(nucleotides),
                'total': sum(nucleotides.values())
            },
            'amino_acids': {
                'count': len(amino_acids),
                'total': sum(amino_acids.values())
            },
            'lipids': {
                'count': len(lipids),
                'total': sum(lipids.values())
            },
        },
        'consumed': consumed,
        'produced': produced
    }
    
    return summary


def export_vbof_to_dict(vbof: VBOFComponents) -> Dict:
    """
    Export VBOF to a dictionary format for saving/loading.
    
    Parameters:
    -----------
    vbof : VBOFComponents
        Complete VBOF
    
    Returns:
    --------
    dict : Serializable dictionary
    """
    return {
        'genome_stoichiometry': vbof.genome_stoichiometry,
        'protein_stoichiometry': vbof.protein_stoichiometry,
        'energy_stoichiometry': vbof.energy_stoichiometry,
        'lipid_stoichiometry': vbof.lipid_stoichiometry,
        'glycan_stoichiometry': vbof.glycan_stoichiometry,
        'combined_stoichiometry': vbof.combined_stoichiometry,
        'metadata': vbof.metadata
    }

