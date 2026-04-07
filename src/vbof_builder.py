"""
VBOF Builder Module
===================

This module assembles all components into a complete Viral Biomass Objective Function for
HMPV and integrates it with the host metabolic model.

Key Functions:
- build_vbof: Assemble complete VBOF stoichiometry
- get_vbof_summary: Summarize VBOF components
- export_vbof_to_dict: Serialize VBOF for saving/loading







Author: Syed Mushahid Hussain
"""

import logging
from typing import Dict, Optional
from dataclasses import dataclass

from .exceptions import VBOFConstructionError
from .glycan_analyzer import calculate_glycan_stoichiometry, calculate_glycosylation_atp
from .energy_calculator import calculate_energy_requirements
from .lipid_calculator import calculate_lipid_stoichiometry

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


def build_vbof(
    genome_stoichiometry: Dict[str, float],
    protein_stoichiometry: Dict[str, float],
    genome_length: int,
    total_amino_acids: int,
    num_proteins: int = 9,
    virion_diameter_nm: float = 209.0,
    include_lipids: bool = True,
    include_energy: bool = True,
    include_glycans: bool = True,
    f_copy_number: int = 350,
    g_copy_number: int = 250,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
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
    
    
    """
    logger.info("Building VBOF...")
    
    # Validate inputs
    if not genome_stoichiometry:
        raise VBOFConstructionError("Genome stoichiometry is empty")
    if not protein_stoichiometry:
        raise VBOFConstructionError("Protein stoichiometry is empty")
    
    glycosylation_atp = (
        calculate_glycosylation_atp(f_copy_number, g_copy_number) if include_glycans else 0.0
    )
    
    # Calculate energy requirements
    if include_energy:
        energy_stoichiometry = calculate_energy_requirements(
            genome_length=genome_length,
            total_amino_acids=total_amino_acids,
            num_proteins=num_proteins,
            glycosylation_atp=glycosylation_atp,
        )
    else:
        energy_stoichiometry = {}
        logger.info("Energy stoichiometry skipped (include_energy=False)")
    
    # Calculate lipid requirements (morphology-aware surface area)
    lipid_stoichiometry = calculate_lipid_stoichiometry(
        virion_diameter_nm=virion_diameter_nm,
        include_lipids=include_lipids,
        morphology=morphology,
        virion_length_nm=virion_length_nm,
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
        'virion_length_nm': virion_length_nm,
        'morphology': morphology,
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

