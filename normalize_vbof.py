#!/usr/bin/env python3
"""
VBOF Normalization Script
=========================

This script normalizes the HMPV VBOF coefficients to mmol per gram dry cell weight (gDCW)
of virion production, which is the standard unit for flux balance analysis.

The raw VBOF has coefficients in molecule counts which are too large for meaningful FBA.

Normalization Formula:
    coefficient_normalized (mmol/gDCW) = coefficient_raw / (total_mass_g * 1000 / MW)

Where:
    - total_mass_g: total mass of one virion (in grams)
    - MW: molecular weight of each component

Usage:
------
    python normalize_vbof.py

Author: Syed Mushahid Hussain
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.config import (
    VBOF_JSON_PATH,
    VBOF_NORMALIZED_JSON_PATH
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Molecular weights (g/mol) - standard biochemistry values
# Source: NIST Chemistry WebBook, Biochemistry textbooks
MOLECULAR_WEIGHTS = {
    # Nucleotides (NTPs)
    'atp_c': 507.18,   # ATP
    'gtp_c': 523.18,   # GTP
    'ctp_c': 483.16,   # CTP
    'utp_c': 484.14,   # UTP
    
    # Amino acids (average ~110 g/mol, using exact values)
    'ala__L_c': 89.09,
    'arg__L_c': 174.20,
    'asn__L_c': 132.12,
    'asp__L_c': 133.10,
    'cys__L_c': 121.16,
    'glu__L_c': 147.13,
    'gln__L_c': 146.15,
    'gly_c': 75.07,
    'his__L_c': 155.16,
    'ile__L_c': 131.17,
    'leu__L_c': 131.17,
    'lys__L_c': 146.19,
    'met__L_c': 149.21,
    'phe__L_c': 165.19,
    'pro__L_c': 115.13,
    'ser__L_c': 105.09,
    'thr__L_c': 119.12,
    'trp__L_c': 204.23,
    'tyr__L_c': 181.19,
    'val__L_c': 117.15,
    
    # Energy carriers
    'adp_c': 427.20,
    'amp_c': 347.22,
    'gdp_c': 443.20,
    'pi_c': 95.98,     # Phosphate
    'ppi_c': 174.95,   # Pyrophosphate
    'h2o_c': 18.02,
    'h_c': 1.01,
    
    # Lipids (approximate average)
    'pc_hs_c': 760.0,      # Phosphatidylcholine
    'pe_hs_c': 720.0,      # Phosphatidylethanolamine
    'ps_hs_c': 780.0,      # Phosphatidylserine
    'sphmyln_hs_c': 730.0, # Sphingomyelin
    'chsterol_c': 386.65,  # Cholesterol
    
    # Glycan precursors
    'uacgam_c': 607.35,    # UDP-N-acetylglucosamine
    'gdpmann_c': 605.34,   # GDP-mannose
    'udpgal_c': 566.30,    # UDP-galactose
    'cmpacna_c': 614.39,   # CMP-N-acetylneuraminic acid
    'gdpfuc_c': 589.33,    # GDP-fucose
    'udpgalfur_c': 607.35, # UDP-N-acetylgalactosamine
}

# Avogadro's number
AVOGADRO = 6.022e23


def calculate_virion_mass(vbof_stoichiometry: dict) -> float:
    """
    Calculate the total mass of one virion based on VBOF stoichiometry.
    
    Parameters:
    -----------
    vbof_stoichiometry : dict
        Dictionary of metabolite IDs to coefficients (raw counts)
    
    Returns:
    --------
    float : Mass of one virion in grams
    """
    total_mass = 0.0
    
    for met_id, coef in vbof_stoichiometry.items():
        if coef < 0:  # Only count consumed metabolites (substrates)
            mw = MOLECULAR_WEIGHTS.get(met_id, 100.0)  # Default MW if not found
            
            # Mass = (molecules) * (MW g/mol) / Avogadro
            mass = abs(coef) * mw / AVOGADRO
            total_mass += mass
    
    return total_mass


def normalize_vbof(vbof_stoichiometry: dict, target_mass_g: float = 1.0) -> dict:
    """
    Normalize VBOF coefficients to mmol per gram dry cell weight.
    
    The standard approach is to normalize such that producing 1 gDCW of 
    virion biomass requires the specified amounts of substrates.
    
    Parameters:
    -----------
    vbof_stoichiometry : dict
        Dictionary of metabolite IDs to raw coefficients
    target_mass_g : float
        Target mass in grams (default: 1.0 gDCW)
    
    Returns:
    --------
    dict : Normalized stoichiometry (mmol/gDCW)
    """
    # Calculate mass of one "unit" based on current coefficients
    virion_mass = calculate_virion_mass(vbof_stoichiometry)
    
    logger.info(f"Calculated virion mass: {virion_mass:.4e} g")
    
    # Calculate scaling factor
    # We need to scale from "per virion" to "per gDCW"
    # Number of virions in 1 gDCW = 1 / virion_mass
    virions_per_gdcw = target_mass_g / virion_mass if virion_mass > 0 else 1.0
    
    logger.info(f"Virions per gDCW: {virions_per_gdcw:.4e}")
    
    normalized = {}
    for met_id, coef in vbof_stoichiometry.items():
        mw = MOLECULAR_WEIGHTS.get(met_id, 100.0)
        
        # Convert from molecules to mmol
        # mmol = molecules / Avogadro * 1000
        # Then scale to per gDCW
        
        # Molecules per gDCW
        molecules_per_gdcw = abs(coef) * virions_per_gdcw
        
        # Convert to mmol/gDCW
        mmol_per_gdcw = molecules_per_gdcw / AVOGADRO * 1000
        
        # Keep the sign
        if coef < 0:
            normalized[met_id] = -mmol_per_gdcw
        else:
            normalized[met_id] = mmol_per_gdcw
    
    return normalized


def normalize_simple(vbof_stoichiometry: dict, normalization_factor: float = None) -> dict:
    """
    Simple normalization: divide all coefficients by a common factor.
    
    This approach normalizes by the largest coefficient to get values in 
    a reasonable range for FBA (typically -1 to +1 for major components).
    
    Parameters:
    -----------
    vbof_stoichiometry : dict
        Dictionary of metabolite IDs to raw coefficients
    normalization_factor : float
        Factor to divide by. If None, uses max absolute coefficient.
    
    Returns:
    --------
    dict : Normalized stoichiometry
    """
    if normalization_factor is None:
        # Find the largest coefficient (by absolute value)
        normalization_factor = max(abs(c) for c in vbof_stoichiometry.values())
    
    logger.info(f"Normalization factor: {normalization_factor:.4e}")
    
    normalized = {
        met_id: coef / normalization_factor 
        for met_id, coef in vbof_stoichiometry.items()
    }
    
    return normalized


def main():
    """
    Main function to normalize the VBOF.
    """
    logger.info("=" * 70)
    logger.info("VBOF NORMALIZATION")
    logger.info("=" * 70)
    
    # Load raw VBOF
    logger.info(f"\nLoading VBOF from: {VBOF_JSON_PATH}")
    with open(VBOF_JSON_PATH, 'r') as f:
        vbof_data = json.load(f)
    
    raw_stoichiometry = vbof_data['combined_stoichiometry']
    
    logger.info(f"Raw coefficients (sample):")
    for i, (k, v) in enumerate(list(raw_stoichiometry.items())[:5]):
        logger.info(f"  {k}: {v:.4e}")
    logger.info(f"  ...")
    
    # Calculate total consumption
    total_consumed = sum(abs(c) for c in raw_stoichiometry.values() if c < 0)
    logger.info(f"\nTotal raw consumed: {total_consumed:.4e}")
    
    # =========================================================================
    # Normalization Method 1: Simple (divide by total consumed)
    # =========================================================================
    logger.info("\n" + "-" * 50)
    logger.info("NORMALIZATION METHOD: Simple (per unit flux)")
    logger.info("-" * 50)
    
    # Use total consumed as normalization factor
    # This gives coefficients that represent fraction of total
    normalized = normalize_simple(raw_stoichiometry, normalization_factor=total_consumed)
    
    logger.info(f"\nNormalized coefficients (sample):")
    for i, (k, v) in enumerate(list(normalized.items())[:5]):
        logger.info(f"  {k}: {v:.6f}")
    logger.info(f"  ...")
    
    # Verify normalization
    total_normalized_consumed = sum(abs(c) for c in normalized.values() if c < 0)
    total_normalized_produced = sum(c for c in normalized.values() if c > 0)
    logger.info(f"\nAfter normalization:")
    logger.info(f"  Total consumed: {total_normalized_consumed:.4f}")
    logger.info(f"  Total produced: {total_normalized_produced:.4f}")
    
    # Save normalized VBOF
    normalized_vbof = {
        'combined_stoichiometry': normalized,
        'metadata': {
            **vbof_data['metadata'],
            'normalization_method': 'simple',
            'normalization_factor': total_consumed,
            'normalized_at': datetime.now().isoformat()
        },
        # Keep original component stoichiometries for reference
        'raw_genome_stoichiometry': vbof_data.get('genome_stoichiometry', {}),
        'raw_protein_stoichiometry': vbof_data.get('protein_stoichiometry', {}),
        'raw_energy_stoichiometry': vbof_data.get('energy_stoichiometry', {}),
        'raw_lipid_stoichiometry': vbof_data.get('lipid_stoichiometry', {}),
        'raw_glycan_stoichiometry': vbof_data.get('glycan_stoichiometry', {}),
    }
    
    logger.info(f"\nSaving normalized VBOF to: {VBOF_NORMALIZED_JSON_PATH}")
    with open(VBOF_NORMALIZED_JSON_PATH, 'w') as f:
        json.dump(normalized_vbof, f, indent=2)
    
    logger.info("Done!")
    
    # Print summary
    print("\n" + "=" * 70)
    print("NORMALIZATION SUMMARY")
    print("=" * 70)
    print(f"\nInput:  {VBOF_JSON_PATH}")
    print(f"Output: {VBOF_NORMALIZED_JSON_PATH}")
    print(f"\nNormalization factor: {total_consumed:.4e}")
    print(f"\nKey normalized coefficients:")
    key_mets = ['atp_c', 'gtp_c', 'leu__L_c', 'h2o_c', 'pi_c']
    for met in key_mets:
        if met in normalized:
            print(f"  {met}: {normalized[met]:.6f}")
    
    return normalized


if __name__ == "__main__":
    normalized = main()

