"""
Energy Calculator Module
========================

Calculates ATP and GTP requirements for HMPV viral replication,
including RNA synthesis, protein translation, and glycosylation costs.

Key Functions:
- calculate_energy_requirements: Calculate ATP/GTP costs for replication

Author: Syed Mushahid Hussain
"""

import logging
from typing import Dict

from .exceptions import EnergyDataMissingError

logger = logging.getLogger(__name__)

ENERGY_COSTS = {
    'rna_polymerization_atp_per_nt': 2,  # 2 ATP equivalents per nucleotide
    'trna_charging_atp_per_aa': 2,        # Amino acid activation
    'elongation_gtp_per_aa': 2,           # EF-Tu and EF-G 
    'initiation_gtp': 1,                  # Initiation factor IF-2
    'termination_gtp': 1,                 # Release factor RF-3
}


def calculate_energy_requirements(
    genome_length: int,
    total_amino_acids: int,
    num_proteins: int = 9,
    glycosylation_atp: float = 0.0,
) -> Dict[str, float]:
    """
    Calculate ATP and GTP requirements for viral replication.

    Energy costs include:
    1. RNA synthesis: NTP → NMP + PPi (2 ATP equiv per nt)
    2. tRNA charging: AA + ATP → AA-tRNA + AMP + PPi (2 ATP per aa)
    3. Translation elongation: 2 GTP per amino acid
    4. Initiation/termination: additional GTP costs
    5. Glycoprotein glycosylation: ATP for F/G N- and O-linked transfer (optional)

    Parameters:
    -----------
    genome_length : int
        Length of viral genome in nucleotides
    total_amino_acids : int
        Total amino acids in all proteins × copy numbers
    num_proteins : int
        Number of distinct proteins (for init/term calculation)
    glycosylation_atp : float
        Extra ATP for glycosylation machinery (from calculate_glycosylation_atp); 0 to omit

    Returns:
    --------
    dict : Energy stoichiometry (negative = consumption)

    Note:
    -----
    These are minimum estimates. Actual costs may be higher due to
    proofreading, chaperones, and other processes.

    Raises:
    -------
    EnergyDataMissingError : If input values are invalid
    
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

    # Total energy (glycosylation ATP → ADP + Pi, same bookkeeping as RNA synthesis ATP)
    total_atp = rna_energy_atp + protein_atp + glycosylation_atp
    total_gtp = protein_gtp + init_term_gtp


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
    if glycosylation_atp:
        logger.info(f"  Glycosylation ATP: {glycosylation_atp}")
    logger.info(f"  Translation GTP: {protein_gtp}")
    logger.info(f"  Total ATP: {total_atp}")
    logger.info(f"  Total GTP: {total_gtp}")

    return stoichiometry
