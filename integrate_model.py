#!/usr/bin/env python3
"""
HMPV Model Integration Script
=============================

This script integrates the HMPV VBOF into the host metabolic model.

Supported Host Models:
----------------------
- model_clean.xml: iHBEC model with R_biomass_hbec (default)
- iHsaEC21_clean.xml: Human airway epithelial cell model (alternative)

Usage:
------
    python integrate_model.py

Input:
------
    - output_HBEC/hmpv_vbof_normalized.json: VBOF stoichiometry from build_vbof.py
    - Data/smbl/model_clean.xml: Host metabolic model (iHBEC)

Output:
-------
    - output_HBEC/{model_id}_with_HMPV_VBOF.xml: Integrated model
    - output_HBEC/integration_summary.txt: Integration summary
    - Console output with validation results

Author: Syed Mushahid Hussain
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.model_integration import (
    load_host_model,
    integrate_vbof,
    validate_integrated_model,
    save_integrated_model,
    get_integration_summary,
    map_vbof_metabolites
)
from src.config import (
    MODEL_DIR,
    DEFAULT_HOST_MODEL,
    OUTPUT_DIR,
    VBOF_REACTION_ID,
    VBOF_NORMALIZED_JSON_PATH,
    INTEGRATION_SUMMARY_PATH,
    INTEGRATED_MODEL_XML_SUFFIX
)
from src.exceptions import HMPVModelError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
logger = logging.getLogger(__name__)


def load_vbof_from_json(vbof_path: str) -> dict:
    """
    Load VBOF stoichiometry from JSON file.
    
    Parameters:
    -----------
    vbof_path : str
        Path to VBOF JSON file
    
    Returns:
    --------
    dict : VBOF data including combined_stoichiometry
    """
    vbof_path = Path(vbof_path)
    
    if not vbof_path.exists():
        raise FileNotFoundError(f"VBOF file not found: {vbof_path}")
    
    with open(vbof_path, 'r') as f:
        vbof_data = json.load(f)
    
    logger.info(f"Loaded VBOF from: {vbof_path}")
    logger.info(f"  Total metabolites: {vbof_data['metadata']['total_metabolites']}")
    
    return vbof_data


def main():
    """
    Main integration pipeline.
    
    Steps:
    1. Load VBOF stoichiometry from JSON
    2. Load host metabolic model
    3. Map metabolites
    4. Integrate VBOF
    5. Validate integrated model
    6. Save results
    """
    logger.info("=" * 70)
    logger.info("HMPV MODEL INTEGRATION PIPELINE")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # =====================================================================
        # STEP 1: Load VBOF from JSON
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 1: LOADING VBOF")
        logger.info("=" * 70)
        
        vbof_data = load_vbof_from_json(VBOF_NORMALIZED_JSON_PATH)
        vbof_stoichiometry = vbof_data['combined_stoichiometry']
        
        logger.info(f"\nVBOF Components:")
        logger.info(f"  Consumed metabolites: {sum(1 for v in vbof_stoichiometry.values() if v < 0)}")
        logger.info(f"  Produced metabolites: {sum(1 for v in vbof_stoichiometry.values() if v > 0)}")
        
        # =====================================================================
        # STEP 2: Load Host Model
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 2: LOADING HOST MODEL")
        logger.info("=" * 70)
        
        model_path = MODEL_DIR / DEFAULT_HOST_MODEL
        host_model = load_host_model(model_path)
        
        logger.info(f"\nHost Model Summary:")
        logger.info(f"  Model ID: {host_model.id}")
        logger.info(f"  Reactions: {len(host_model.reactions)}")
        logger.info(f"  Metabolites: {len(host_model.metabolites)}")
        logger.info(f"  Genes: {len(host_model.genes)}")
        
        # =====================================================================
        # STEP 3: Preview Metabolite Mapping
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 3: METABOLITE MAPPING PREVIEW")
        logger.info("=" * 70)
        
        mappings, unmapped = map_vbof_metabolites(host_model, vbof_stoichiometry)
        
        found_count = sum(1 for m in mappings.values() if m.found)
        logger.info(f"\nMapping Results:")
        logger.info(f"  Found: {found_count}/{len(mappings)} metabolites")
        logger.info(f"  Unmapped: {len(unmapped)}")
        
        if unmapped:
            logger.warning(f"\nUnmapped metabolites (will be skipped):")
            for uid in unmapped:
                m = mappings[uid]
                logger.warning(f"  {uid}")
                if m.suggestions:
                    logger.warning(f"    Suggestions: {m.suggestions[:3]}")
        
        # =====================================================================
        # STEP 4: Integrate VBOF
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 4: INTEGRATING VBOF INTO HOST MODEL")
        logger.info("=" * 70)
        
        integrated_model, vbof_reaction, unmapped_final = integrate_vbof(
            model=host_model,
            vbof_stoichiometry=vbof_stoichiometry,
            reaction_id=VBOF_REACTION_ID,
            skip_unmapped=True
        )
        
        logger.info(f"\nIntegration Complete:")
        logger.info(f"  VBOF Reaction ID: {vbof_reaction.id}")
        logger.info(f"  Metabolites in VBOF: {len(vbof_reaction.metabolites)}")
        logger.info(f"  Skipped (unmapped): {len(unmapped_final)}")
        
        # =====================================================================
        # STEP 5: Validate Integrated Model
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 5: VALIDATING INTEGRATED MODEL")
        logger.info("=" * 70)
        
        validation = validate_integrated_model(integrated_model, VBOF_REACTION_ID)
        
        logger.info(f"\nValidation Results:")
        logger.info(f"  Valid: {validation['valid']}")
        
        if validation['fba_result']:
            logger.info(f"  FBA Status: {validation['fba_result']['status']}")
            if validation['fba_result']['objective_value']:
                logger.info(f"  Max VBOF Flux: {validation['fba_result']['objective_value']:.6f}")
        
        if validation['errors']:
            for error in validation['errors']:
                logger.error(f"  ERROR: {error}")
        
        if validation['warnings']:
            logger.info(f"\n  Warnings ({len(validation['warnings'])}):")
            for warning in validation['warnings'][:10]:  # Limit output
                logger.warning(f"    {warning}")
            if len(validation['warnings']) > 10:
                logger.warning(f"    ... and {len(validation['warnings']) - 10} more")
        
        # =====================================================================
        # STEP 6: Save Results
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 6: SAVING RESULTS")
        logger.info("=" * 70)
        
        # Save integrated model (SBML)
        model_output_path = save_integrated_model(
            model=integrated_model,
            output_path=OUTPUT_DIR / f"{integrated_model.id}{INTEGRATED_MODEL_XML_SUFFIX}",
            format='sbml'
        )
        logger.info(f"Saved integrated model: {model_output_path}")
        
        # Generate and save integration summary
        summary = get_integration_summary(
            model=integrated_model,
            vbof_reaction_id=VBOF_REACTION_ID,
            unmapped=unmapped_final,
            validation=validation
        )
        
        with open(INTEGRATION_SUMMARY_PATH, 'w') as f:
            f.write(summary)
            f.write(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        logger.info(f"Saved integration summary: {INTEGRATION_SUMMARY_PATH}")
        
        # Print summary
        print(summary)
        
        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("INTEGRATION COMPLETE")
        logger.info("=" * 70)
        
        logger.info(f"\nHMPV VBOF successfully integrated into {integrated_model.id}!")
        logger.info(f"\nOutput files:")
        logger.info(f"  1. {model_output_path}")
        logger.info(f"  2. {INTEGRATION_SUMMARY_PATH}")
        
      
        
        return integrated_model, validation
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        logger.error("Please run build_vbof.py first to generate the VBOF")
        sys.exit(1)
    except HMPVModelError as e:
        logger.error(f"Integration failed: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    model, validation = main()

