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
    VBOF_NORMALIZED_JSON_PATH,
    METABOLITE_MOLECULAR_WEIGHTS,
    STRUCTURAL_METABOLITES,
    NUCLEOTIDE_BIGG_IDS,
    AMINO_ACID_BIGG_IDS,
    LIPID_FRACTIONS,
    GLYCAN_PRECURSOR_BIGG_IDS,
    COMMON_METABOLITE_IDS,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Molecular weights and structural metabolites are imported from src.config above.
# Edit those centralized definitions; this file uses them automatically.
# Alias for local readability
MOLECULAR_WEIGHTS = METABOLITE_MOLECULAR_WEIGHTS

# Avogadro's number
AVOGADRO = 6.022e23
# Add genome-only ATP/GTP contribution manually
GENOME_ATP_MOLECULES = 5103  # 5103 A nucleotides in genome
GENOME_GTP_MOLECULES = 2538

def calculate_virion_mass(vbof_stoichiometry: dict) -> float:
    """
    Calculate virion dry mass (grams) from structural metabolites only.

    Only metabolites in STRUCTURAL_METABOLITES contribute.  atp_c / gtp_c are
    handled separately using their genome-nucleotide counts so that the large
    energy-cost coefficients in the VBOF are excluded.
    """
    total_mass = 0.0
    contributions = {}

    for met_id, coef in vbof_stoichiometry.items():
        if met_id in STRUCTURAL_METABOLITES and coef < 0:
            if met_id not in MOLECULAR_WEIGHTS:
                raise KeyError(
                    f"Structural metabolite '{met_id}' has no entry in "
                    f"MOLECULAR_WEIGHTS — add its MW before running."
                )
            mw = MOLECULAR_WEIGHTS[met_id]
            mass = abs(coef) * mw / AVOGADRO
            contributions[met_id] = mass
            total_mass += mass

    # Genome-only ATP (A nucleotides) and GTP (G nucleotides)
    atp_id = NUCLEOTIDE_BIGG_IDS['A']
    gtp_id = NUCLEOTIDE_BIGG_IDS['G']
    atp_mass = GENOME_ATP_MOLECULES * MOLECULAR_WEIGHTS[atp_id] / AVOGADRO
    gtp_mass = GENOME_GTP_MOLECULES * MOLECULAR_WEIGHTS[gtp_id] / AVOGADRO
    contributions[f'{atp_id} (genome only)'] = atp_mass
    contributions[f'{gtp_id} (genome only)'] = gtp_mass
    total_mass += atp_mass + gtp_mass

    # Log every contribution for verification
    logger.info("  Metabolite mass contributions:")
    for met, mass in sorted(contributions.items(), key=lambda x: -x[1]):
        logger.info(f"    {met:<30s}: {mass:.4e} g  ({mass / total_mass * 100:.2f}%)")

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

    # Calculate and log virion mass
    virion_mass = calculate_virion_mass(raw_stoichiometry)
    virions_per_gdcw = 1.0 / virion_mass if virion_mass > 0 else 0
    logger.info(f"Calculated virion mass: {virion_mass:.4e} g")
    logger.info(f"Virions per gDCW:       {virions_per_gdcw:.4e}")
    
    # =========================================================================
    # Normalization Method 1: Simple (divide by total consumed)
    # =========================================================================
    logger.info("\n" + "-" * 50)
    logger.info("NORMALIZATION METHOD: Simple (per unit flux)")
    logger.info("-" * 50)
    
    # Use total consumed as normalization factor
    # This gives coefficients that represent fraction of total
    #normalized = normalize_simple(raw_stoichiometry, normalization_factor=total_consumed)
    normalized = normalize_vbof(raw_stoichiometry, target_mass_g=1.0)
    
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
    key_mets = [
        NUCLEOTIDE_BIGG_IDS['A'],
        NUCLEOTIDE_BIGG_IDS['G'],
        AMINO_ACID_BIGG_IDS['L'],
        COMMON_METABOLITE_IDS['h2o'],
        COMMON_METABOLITE_IDS['pi'],
    ]
    for met in key_mets:
        if met in normalized:
            print(f"  {met}: {normalized[met]:.6f}")
    
    # =========================================================================
    # Verification: corrected virion mass and category breakdown
    # =========================================================================
    print("\n" + "=" * 70)
    print("VIRION MASS VERIFICATION")
    print("=" * 70)
    print(f"\nCorrected virion mass: {virion_mass:.4e} g  "
          f"({virion_mass * 1e15:.2f} femtograms)")

    # --- per-metabolite contributions (mirroring calculate_virion_mass logic) ---
    atp_id = NUCLEOTIDE_BIGG_IDS['A']
    gtp_id = NUCLEOTIDE_BIGG_IDS['G']
    contributions = {}
    for met_id, coef in raw_stoichiometry.items():
        if met_id in STRUCTURAL_METABOLITES and coef < 0:
            mw = MOLECULAR_WEIGHTS[met_id]
            contributions[met_id] = abs(coef) * mw / AVOGADRO
    contributions[f'{atp_id} (genome only)'] = GENOME_ATP_MOLECULES * MOLECULAR_WEIGHTS[atp_id] / AVOGADRO
    contributions[f'{gtp_id} (genome only)'] = GENOME_GTP_MOLECULES * MOLECULAR_WEIGHTS[gtp_id] / AVOGADRO

    # Top 10 by mass
    top10 = sorted(contributions.items(), key=lambda x: -x[1])[:10]
    print("\nTop 10 metabolites by mass contribution:")
    for met, mass in top10:
        print(f"  {met:<30s}: {mass:.4e} g  ({mass / virion_mass * 100:.1f}%)")

    # Category breakdown — IDs sourced from config variables
    nucleotide_ids = {
        NUCLEOTIDE_BIGG_IDS['C'],
        NUCLEOTIDE_BIGG_IDS['U'],
        f'{atp_id} (genome only)',
        f'{gtp_id} (genome only)',
    }
    categories = {
        'Amino acids': set(AMINO_ACID_BIGG_IDS.values()),
        'Lipids':      set(LIPID_FRACTIONS.keys()),
        'Glycans':     set(GLYCAN_PRECURSOR_BIGG_IDS.values()),
        'Nucleotides': nucleotide_ids,
    }
    print("\nMass breakdown by category:")
    for cat, ids in categories.items():
        cat_mass = sum(contributions.get(i, 0.0) for i in ids)
        print(f"  {cat:<14s}: {cat_mass:.4e} g  ({cat_mass / virion_mass * 100:.1f}%)")

    return normalized


if __name__ == "__main__":
    normalized = main()

