#!/usr/bin/env python3
"""
Clean Host Model Script
=======================

This script removes the SARS-CoV-2 VBOF from the iHsaEC21 model to create
a clean host-only model for HMPV integration.

The original model (iHsaEC21_PLUS_SARS_CoV_2) contains a pre-integrated
SARS-CoV-2 VBOF which would interfere with HMPV-specific analysis.

Usage:
------
    python clean_host_model.py

Input:
------
    - Data/smbl/iHsaEC21.xml (contains SARS-CoV-2 VBOF)

Output:
-------
    - Data/smbl/iHsaEC21_clean.xml (pure host model)
    - output/model_cleaning_report.txt

Author: Syed Mushahid Hussain
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

import cobra
from cobra.io import read_sbml_model, write_sbml_model

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.config import (
    MODEL_DIR,
    HOST_MODEL_ORIGINAL_PATH,
    HOST_MODEL_CLEAN_PATH,
    MODEL_CLEANING_REPORT_PATH
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def find_viral_reactions(model):
    """
    Find all viral-related reactions in the model.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    
    Returns:
    --------
    list : List of viral reaction IDs
    """
    viral_keywords = [
        'vbof', 'viral', 'sars', 'cov', 'covid', 
        'corona', 'virus', 'virion'
    ]
    
    viral_reactions = []
    for reaction in model.reactions:
        rid_lower = reaction.id.lower()
        rname_lower = (reaction.name or '').lower()
        
        for keyword in viral_keywords:
            if keyword in rid_lower or keyword in rname_lower:
                viral_reactions.append(reaction.id)
                break
    
    return viral_reactions


def analyze_model(model):
    """
    Analyze the model and return summary information.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    
    Returns:
    --------
    dict : Model summary
    """
    summary = {
        'model_id': model.id,
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'objective': str(model.objective.expression),
        'viral_reactions': find_viral_reactions(model)
    }
    return summary


def remove_sars_cov2_vbof(model, vbof_id='VBOF'):
    """
    Remove the SARS-CoV-2 VBOF from the model.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model with SARS-CoV-2 VBOF
    vbof_id : str
        The reaction ID of the SARS-CoV-2 VBOF
    
    Returns:
    --------
    cobra.Model : Cleaned model
    dict : Details of what was removed
    """
    removed_info = {
        'reaction_id': None,
        'reaction_name': None,
        'metabolites': [],
        'success': False
    }
    
    # Check if VBOF exists
    if vbof_id not in [r.id for r in model.reactions]:
        logger.warning(f"VBOF reaction '{vbof_id}' not found in model")
        return model, removed_info
    
    # Get VBOF details before removal
    vbof = model.reactions.get_by_id(vbof_id)
    removed_info['reaction_id'] = vbof.id
    removed_info['reaction_name'] = vbof.name
    removed_info['metabolites'] = [(m.id, c) for m, c in vbof.metabolites.items()]
    removed_info['bounds'] = (vbof.lower_bound, vbof.upper_bound)
    
    logger.info(f"Removing SARS-CoV-2 VBOF: {vbof_id}")
    logger.info(f"  Metabolites in VBOF: {len(vbof.metabolites)}")
    
    # Remove the reaction
    model.remove_reactions([vbof_id])
    
    # Update model ID to reflect it's now clean
    model.id = model.id.replace('_PLUS_SARS_CoV_2', '_CLEAN')
    if '_CLEAN' not in model.id:
        model.id = model.id + '_CLEAN'
    
    removed_info['success'] = True
    logger.info(f"VBOF removed successfully!")
    logger.info(f"New model ID: {model.id}")
    
    return model, removed_info


def set_default_objective(model):
    """
    Set a default objective for the cleaned model.
    
    For a host cell model without viral components, we typically
    don't set an objective until we add our own VBOF.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    """
    # Check if there's a biomass reaction we can use
    biomass_keywords = ['biomass', 'growth', 'maintenance']
    
    for reaction in model.reactions:
        rid_lower = reaction.id.lower()
        for keyword in biomass_keywords:
            if keyword in rid_lower:
                logger.info(f"Found potential objective: {reaction.id}")
                # Don't automatically set it, just report
                return
    
    logger.info("No default biomass objective found - model is ready for VBOF integration")


def generate_report(original_summary, cleaned_summary, removed_info, output_path):
    """
    Generate a report of the cleaning process.
    
    Parameters:
    -----------
    original_summary : dict
        Summary of original model
    cleaned_summary : dict
        Summary of cleaned model
    removed_info : dict
        Information about what was removed
    output_path : str
        Path to save the report
    """
    report = f"""
================================================================================
HOST MODEL CLEANING REPORT
================================================================================

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

ORIGINAL MODEL:
---------------
Model ID: {original_summary['model_id']}
Reactions: {original_summary['reactions']}
Metabolites: {original_summary['metabolites']}
Genes: {original_summary['genes']}
Objective: {original_summary['objective']}

Viral Reactions Found:
"""
    for vr in original_summary['viral_reactions']:
        report += f"  - {vr}\n"
    
    report += f"""
REMOVED SARS-CoV-2 VBOF:
------------------------
Reaction ID: {removed_info['reaction_id']}
Reaction Name: {removed_info['reaction_name']}
Bounds: {removed_info.get('bounds', 'N/A')}
Number of metabolites: {len(removed_info['metabolites'])}

Metabolites in SARS-CoV-2 VBOF:
"""
    for mid, coef in removed_info['metabolites']:
        report += f"  {mid}: {coef}\n"
    
    report += f"""
CLEANED MODEL:
--------------
Model ID: {cleaned_summary['model_id']}
Reactions: {cleaned_summary['reactions']}
Metabolites: {cleaned_summary['metabolites']}
Genes: {cleaned_summary['genes']}

Viral Reactions After Cleaning:
"""
    if cleaned_summary['viral_reactions']:
        for vr in cleaned_summary['viral_reactions']:
            report += f"  - {vr}\n"
    else:
        report += "  None - Model is clean!\n"
    
    report += f"""
CHANGES:
--------
Reactions removed: {original_summary['reactions'] - cleaned_summary['reactions']}
Model ready for HMPV VBOF integration: YES

================================================================================
"""
    
    with open(output_path, 'w') as f:
        f.write(report)
    
    logger.info(f"Report saved to: {output_path}")
    return report


def main():
    """
    Main function to clean the host model.
    """
    logger.info("=" * 70)
    logger.info("HOST MODEL CLEANING SCRIPT")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Ensure output directory exists
    MODEL_CLEANING_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    
    # =========================================================================
    # STEP 1: Load original model
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 1: LOADING ORIGINAL MODEL")
    logger.info("=" * 70)
    
    if not HOST_MODEL_ORIGINAL_PATH.exists():
        logger.error(f"Model file not found: {HOST_MODEL_ORIGINAL_PATH}")
        sys.exit(1)
    
    logger.info(f"Loading model from: {HOST_MODEL_ORIGINAL_PATH}")
    model = read_sbml_model(str(HOST_MODEL_ORIGINAL_PATH))
    
    # Analyze original model
    original_summary = analyze_model(model)
    
    logger.info(f"\nOriginal Model Summary:")
    logger.info(f"  Model ID: {original_summary['model_id']}")
    logger.info(f"  Reactions: {original_summary['reactions']}")
    logger.info(f"  Metabolites: {original_summary['metabolites']}")
    logger.info(f"  Genes: {original_summary['genes']}")
    logger.info(f"  Viral reactions: {len(original_summary['viral_reactions'])}")
    
    for vr in original_summary['viral_reactions']:
        logger.info(f"    - {vr}")
    
    # =========================================================================
    # STEP 2: Remove SARS-CoV-2 VBOF
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: REMOVING SARS-CoV-2 VBOF")
    logger.info("=" * 70)
    
    model, removed_info = remove_sars_cov2_vbof(model, vbof_id='VBOF')
    
    if not removed_info['success']:
        logger.warning("VBOF was not removed - model may already be clean")
    
    # Set default objective
    set_default_objective(model)
    
    # =========================================================================
    # STEP 3: Validate cleaned model
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: VALIDATING CLEANED MODEL")
    logger.info("=" * 70)
    
    cleaned_summary = analyze_model(model)
    
    logger.info(f"\nCleaned Model Summary:")
    logger.info(f"  Model ID: {cleaned_summary['model_id']}")
    logger.info(f"  Reactions: {cleaned_summary['reactions']}")
    logger.info(f"  Metabolites: {cleaned_summary['metabolites']}")
    logger.info(f"  Genes: {cleaned_summary['genes']}")
    logger.info(f"  Viral reactions: {len(cleaned_summary['viral_reactions'])}")
    
    if cleaned_summary['viral_reactions']:
        logger.warning("Some viral reactions remain:")
        for vr in cleaned_summary['viral_reactions']:
            logger.warning(f"    - {vr}")
    else:
        logger.info("  ✓ No viral reactions - model is clean!")
    
    # =========================================================================
    # STEP 4: Save cleaned model
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 4: SAVING CLEANED MODEL")
    logger.info("=" * 70)
    
    logger.info(f"Saving cleaned model to: {HOST_MODEL_CLEAN_PATH}")
    write_sbml_model(model, str(HOST_MODEL_CLEAN_PATH))
    logger.info("Model saved successfully!")
    
    # Generate report
    report = generate_report(
        original_summary, 
        cleaned_summary, 
        removed_info, 
        MODEL_CLEANING_REPORT_PATH
    )
    
    # Print summary
    print(report)
    
    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("CLEANING COMPLETE")
    logger.info("=" * 70)
    
    logger.info(f"\nOutput files:")
    logger.info(f"  1. {HOST_MODEL_CLEAN_PATH} (clean host model)")
    logger.info(f"  2. {MODEL_CLEANING_REPORT_PATH} (cleaning report)")
    
  
    
    return model


if __name__ == "__main__":
    model = main()

