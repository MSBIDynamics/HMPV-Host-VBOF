#!/usr/bin/env python3
"""
Dual-Objective Gene Knockout Analysis
======================================

This script performs comprehensive gene and reaction knockout analysis comparing host cell and viral growth
for each gene and reaction deletion.


Methodology:
------------
1. Load integrated model (host + HMPV VBOF)
2. Run gene knockouts with host BOF as objective -> host growth results
3. Run gene knockouts with VBOF as objective -> virus growth results
4. Merge results into comparison table with percentage calculations
5. Apply configurable thresholds to identify selective targets
6. Optionally run combined objective analysis


Output Files:
-------------
- host_growth_knockout_results.csv: Gene knockouts with BOF objective
- virus_growth_knockout_results.csv: Gene knockouts with VBOF objective
- merged_knockout_comparison.csv: Combined comparison table
- selective_antiviral_targets.csv: Selective targets 
- critical_viral_targets.csv: Critical targets 
- combined_objective_results.csv: Combined objective analysis 
- dual_objective_report.txt: Comprehensive text report

Usage:
------
    python dual_objective_analysis.py [options]

Options:
--------
    --virus-max FLOAT   Maximum virus growth % for selective targets (default: 50)
    --host-min FLOAT    Minimum host growth % for selective targets (default: 80)
    --critical FLOAT    Maximum virus growth % for critical targets (default: 5)
    --combined          Run combined objective analysis
    --bof-weight FLOAT  BOF weight for combined analysis (default: 0.9)
    --vbof-weight FLOAT VBOF weight for combined analysis (default: 0.1)

Examples:
---------
    # Default analysis (virus <50%, host >80%)
    python dual_objective_analysis.py
    
    # Strict thresholds (virus <10%, host >90%)
    python dual_objective_analysis.py --virus-max 10 --host-min 90
    
    # With combined objective analysis
    python dual_objective_analysis.py --combined --bof-weight 0.9 --vbof-weight 0.1

Author: Syed Mushahid Hussain
"""

import argparse
import logging
import sys
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional
import warnings

import cobra
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.config import (
    OUTPUT_DIR,
    MODEL_DIR,
    DEFAULT_HOST_MODEL,
    DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR,
    HOST_GROWTH_RESULTS_PATH,
    VIRUS_GROWTH_RESULTS_PATH,
    MERGED_KNOCKOUT_RESULTS_PATH,
    SELECTIVE_TARGETS_PATH,
    CRITICAL_VIRAL_TARGETS_PATH,
    COMBINED_OBJECTIVE_RESULTS_PATH,
    DUAL_OBJECTIVE_REPORT_PATH,
    # Reaction knockout paths
    HOST_REACTION_KNOCKOUT_RESULTS_PATH,
    VIRUS_REACTION_KNOCKOUT_RESULTS_PATH,
    MERGED_REACTION_KNOCKOUT_RESULTS_PATH,
    SELECTIVE_REACTION_TARGETS_PATH,
    CRITICAL_REACTION_TARGETS_PATH,
    REACTION_SUBSYSTEM_ESSENTIALITY_PATH,
    VBOF_REACTION_ID,
    HOST_BOF_REACTION_ID,
    DEFAULT_THRESHOLDS,
    INTEGRATED_MODEL_XML_SUFFIX,
    VBOF_NORMALIZED_JSON_PATH
)
from src.dual_objective_knockout import (
    DualObjectiveKnockout,
    DualObjectiveConfig,
    load_integrated_model
)
from src.model_integration import (
    load_host_model,
    integrate_vbof,
    validate_integrated_model,
    save_integrated_model
)
from src.exceptions import HMPVModelError

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Dual-objective gene knockout analysis for antiviral target identification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dual_objective_analysis.py                          # Default analysis
  python dual_objective_analysis.py --virus-max 10 --host-min 90  # Strict thresholds
  python dual_objective_analysis.py --combined               # With combined objective
        """
    )
    
    # Threshold options
    parser.add_argument(
        '--virus-max', type=float, default=50.0,
        help='Maximum virus growth %% for selective targets (default: 50)'
    )
    parser.add_argument(
        '--host-min', type=float, default=80.0,
        help='Minimum host growth %% for selective targets (default: 80)'
    )
    parser.add_argument(
        '--critical', type=float, default=5.0,
        help='Maximum virus growth %% for critical targets (default: 5)'
    )
    
    # Combined objective options
    parser.add_argument(
        '--combined', action='store_true',
        help='Run combined objective analysis'
    )
    parser.add_argument(
        '--bof-weight', type=float, default=0.9,
        help='BOF weight for combined analysis (default: 0.9)'
    )
    parser.add_argument(
        '--vbof-weight', type=float, default=0.1,
        help='VBOF weight for combined analysis (default: 0.1)'
    )
    
    # Model options
    parser.add_argument(
        '--model', type=str, default=None,
        help='Path to integrated model (auto-detected if not specified)'
    )
    
    return parser.parse_args()


# =============================================================================
# MODEL LOADING AND INTEGRATION
# =============================================================================

def find_or_create_integrated_model() -> Path:
    """
    Find existing integrated model or create one.
    
    Returns:
    --------
    Path : Path to integrated model file
    """
    # Look for existing integrated model
    integrated_models = list(OUTPUT_DIR.glob(f"*{INTEGRATED_MODEL_XML_SUFFIX}"))
    
    if integrated_models:
        # Use most recent
        model_path = sorted(integrated_models, key=lambda p: p.stat().st_mtime, reverse=True)[0]
        logger.info(f"Found existing integrated model: {model_path.name}")
        return model_path
    
    # Need to create integrated model
    logger.info("No integrated model found. Creating one...")
    
    # Check for VBOF JSON
    if not VBOF_NORMALIZED_JSON_PATH.exists():
        raise FileNotFoundError(
            f"VBOF file not found: {VBOF_NORMALIZED_JSON_PATH}\n"
            "Please run build_vbof.py first to generate the VBOF"
        )
    
    # Load VBOF
    with open(VBOF_NORMALIZED_JSON_PATH, 'r') as f:
        vbof_data = json.load(f)
    vbof_stoichiometry = vbof_data['combined_stoichiometry']
    
    # Load host model
    model_path = MODEL_DIR / DEFAULT_HOST_MODEL
    host_model = load_host_model(model_path)
    
    # Integrate VBOF
    integrated_model, vbof_reaction, unmapped = integrate_vbof(
        model=host_model,
        vbof_stoichiometry=vbof_stoichiometry,
        reaction_id=VBOF_REACTION_ID,
        skip_unmapped=True
    )
    
    # Save integrated model
    output_path = OUTPUT_DIR / f"{integrated_model.id}{INTEGRATED_MODEL_XML_SUFFIX}"
    save_integrated_model(integrated_model, output_path, format='sbml')
    
    logger.info(f"Created integrated model: {output_path}")
    
    return output_path


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """Main pipeline for dual-objective gene knockout analysis."""
    
    # Parse arguments
    args = parse_arguments()
    
    logger.info("=" * 70)
    logger.info("DUAL-OBJECTIVE GENE KNOCKOUT ANALYSIS")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    try:
        # =====================================================================
        # STEP 1: Load or create integrated model
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 1: LOADING INTEGRATED MODEL")
        logger.info("=" * 70)
        
        if args.model:
            model_path = Path(args.model)
        else:
            model_path = find_or_create_integrated_model()
        
        model = load_integrated_model(model_path)
        
        # =====================================================================
        # STEP 2: Configure analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 2: CONFIGURING ANALYSIS")
        logger.info("=" * 70)
        
        # Update thresholds from command line arguments
        thresholds = DEFAULT_THRESHOLDS.copy()
        thresholds['selective_target'] = {
            'virus_max': args.virus_max / 100.0,
            'host_min': args.host_min / 100.0
        }
        thresholds['lethal_virus'] = {
            'virus_max': args.critical / 100.0
        }
        
        config = DualObjectiveConfig(
            vbof_id=VBOF_REACTION_ID,
            bof_id=HOST_BOF_REACTION_ID,
            thresholds=thresholds
        )
        
        logger.info(f"Thresholds configured:")
        logger.info(f"  Selective target: virus < {args.virus_max}%, host > {args.host_min}%")
        logger.info(f"  Critical target: virus < {args.critical}%")
        
        # =====================================================================
        # STEP 3: Initialize analyzer
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 3: INITIALIZING DUAL-OBJECTIVE ANALYZER")
        logger.info("=" * 70)
        
        analyzer = DualObjectiveKnockout(model, config)
        
        # =====================================================================
        # STEP 4: Run host knockout analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 4: HOST KNOCKOUT ANALYSIS")
        logger.info("=" * 70)
        
        host_results = analyzer.run_host_knockout_analysis()
        host_results.to_csv(HOST_GROWTH_RESULTS_PATH, index=False)
        logger.info(f"Saved: {HOST_GROWTH_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 5: Run virus knockout analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 5: VIRUS KNOCKOUT ANALYSIS")
        logger.info("=" * 70)
        
        virus_results = analyzer.run_virus_knockout_analysis()
        virus_results.to_csv(VIRUS_GROWTH_RESULTS_PATH, index=False)
        logger.info(f"Saved: {VIRUS_GROWTH_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 6: Merge results
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 6: MERGING RESULTS")
        logger.info("=" * 70)
        
        merged_results = analyzer.merge_results(host_results, virus_results)
        merged_results.to_csv(MERGED_KNOCKOUT_RESULTS_PATH, index=False)
        logger.info(f"Saved: {MERGED_KNOCKOUT_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 7: Filter selective targets
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 7: FILTERING SELECTIVE TARGETS")
        logger.info("=" * 70)
        
        selective_targets = analyzer.filter_selective_targets(
            merged_results,
            virus_max=args.virus_max,
            host_min=args.host_min
        )
        selective_targets.to_csv(SELECTIVE_TARGETS_PATH, index=False)
        logger.info(f"Saved: {SELECTIVE_TARGETS_PATH}")
        
        # =====================================================================
        # STEP 8: Filter critical viral targets
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 8: FILTERING CRITICAL VIRAL TARGETS")
        logger.info("=" * 70)
        
        critical_targets = analyzer.filter_critical_viral_targets(
            merged_results,
            virus_max=args.critical
        )
        critical_targets.to_csv(CRITICAL_VIRAL_TARGETS_PATH, index=False)
        logger.info(f"Saved: {CRITICAL_VIRAL_TARGETS_PATH}")
        
        # =====================================================================
        # STEP 9: Combined objective analysis (optional)
        # =====================================================================
        combined_results = None
        if args.combined:
            logger.info("\n" + "=" * 70)
            logger.info("STEP 9: COMBINED OBJECTIVE ANALYSIS")
            logger.info("=" * 70)
            
            combined_results = analyzer.analyze_combined_objective(
                bof_weight=args.bof_weight,
                vbof_weight=args.vbof_weight
            )
            combined_results.to_csv(COMBINED_OBJECTIVE_RESULTS_PATH, index=False)
            logger.info(f"Saved: {COMBINED_OBJECTIVE_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 10: Host REACTION knockout analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 10: HOST REACTION KNOCKOUT ANALYSIS")
        logger.info("=" * 70)
        
        host_rxn_results = analyzer.run_host_reaction_knockout_analysis()
        host_rxn_results.to_csv(HOST_REACTION_KNOCKOUT_RESULTS_PATH, index=False)
        logger.info(f"Saved: {HOST_REACTION_KNOCKOUT_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 11: Virus REACTION knockout analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 11: VIRUS REACTION KNOCKOUT ANALYSIS")
        logger.info("=" * 70)
        
        virus_rxn_results = analyzer.run_virus_reaction_knockout_analysis()
        virus_rxn_results.to_csv(VIRUS_REACTION_KNOCKOUT_RESULTS_PATH, index=False)
        logger.info(f"Saved: {VIRUS_REACTION_KNOCKOUT_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 12: Merge reaction results
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 12: MERGING REACTION RESULTS")
        logger.info("=" * 70)
        
        merged_rxn_results = analyzer.merge_reaction_results(host_rxn_results, virus_rxn_results)
        merged_rxn_results.to_csv(MERGED_REACTION_KNOCKOUT_RESULTS_PATH, index=False)
        logger.info(f"Saved: {MERGED_REACTION_KNOCKOUT_RESULTS_PATH}")
        
        # =====================================================================
        # STEP 13: Filter selective reaction targets
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 13: FILTERING SELECTIVE REACTION TARGETS")
        logger.info("=" * 70)
        
        selective_rxn_targets = analyzer.filter_selective_reaction_targets(
            merged_rxn_results,
            virus_max=args.virus_max,
            host_min=args.host_min
        )
        selective_rxn_targets.to_csv(SELECTIVE_REACTION_TARGETS_PATH, index=False)
        logger.info(f"Saved: {SELECTIVE_REACTION_TARGETS_PATH}")
        
        # =====================================================================
        # STEP 14: Filter critical viral reaction targets
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 14: FILTERING CRITICAL VIRAL REACTION TARGETS")
        logger.info("=" * 70)
        
        critical_rxn_targets = analyzer.filter_critical_reaction_targets(
            merged_rxn_results,
            virus_max=args.critical
        )
        critical_rxn_targets.to_csv(CRITICAL_REACTION_TARGETS_PATH, index=False)
        logger.info(f"Saved: {CRITICAL_REACTION_TARGETS_PATH}")
        
        # =====================================================================
        # STEP 15: Subsystem essentiality analysis
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 15: SUBSYSTEM ESSENTIALITY ANALYSIS")
        logger.info("=" * 70)
        
        subsystem_analysis = analyzer.analyze_subsystem_essentiality(merged_rxn_results)
        subsystem_analysis.to_csv(REACTION_SUBSYSTEM_ESSENTIALITY_PATH, index=False)
        logger.info(f"Saved: {REACTION_SUBSYSTEM_ESSENTIALITY_PATH}")
        
        logger.info(f"Saved: {DUAL_OBJECTIVE_REPORT_PATH}")
        
        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE")
        print("=" * 70)
        
        print(f"\nBaseline Fluxes:")
        print(f"  Host BOF: {analyzer.host_baseline:.6f}")
        print(f"  VBOF: {analyzer.virus_baseline:.6f}")
        
        print(f"\n" + "=" * 50)
        print("GENE KNOCKOUT RESULTS")
        print("=" * 50)
        
        print(f"\n--- SELECTIVE GENE TARGETS ({len(selective_targets)}) ---")
        print(f"(Virus < {args.virus_max}%, Host > {args.host_min}%)")
        for _, row in selective_targets.head(10).iterrows():
            print(f"  {row['gene_id']}: Host {row['host_growth_pct']:.1f}%, Virus {row['virus_growth_pct']:.1f}%")
        
        print(f"\n--- CRITICAL VIRAL GENE TARGETS ({len(critical_targets)}) ---")
        print(f"(Virus < {args.critical}%)")
        for _, row in critical_targets.head(10).iterrows():
            print(f"  {row['gene_id']}: Host {row['host_growth_pct']:.1f}%, Virus {row['virus_growth_pct']:.1f}%")
        
        print(f"\n" + "=" * 50)
        print("REACTION KNOCKOUT RESULTS")
        print("=" * 50)
        
        print(f"\n--- SELECTIVE REACTION TARGETS ({len(selective_rxn_targets)}) ---")
        print(f"(Virus < {args.virus_max}%, Host > {args.host_min}%)")
        for _, row in selective_rxn_targets.head(10).iterrows():
            print(f"  {row['reaction_id']}: Host {row['host_growth_pct']:.1f}%, Virus {row['virus_growth_pct']:.1f}%")
        
        print(f"\n--- CRITICAL VIRAL REACTION TARGETS ({len(critical_rxn_targets)}) ---")
        print(f"(Virus < {args.critical}%)")
        for _, row in critical_rxn_targets.head(10).iterrows():
            print(f"  {row['reaction_id']}: Host {row['host_growth_pct']:.1f}%, Virus {row['virus_growth_pct']:.1f}%")
        
        print(f"\nOutput saved to: {DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR}")
        print("=" * 70)
        
        logger.info("\n" + "=" * 70)
        logger.info("DUAL-OBJECTIVE ANALYSIS COMPLETE")
        logger.info("=" * 70)
        
        return merged_results
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        logger.error("Please ensure all required input files exist")
        sys.exit(1)
    except HMPVModelError as e:
        logger.error(f"Model error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    results = main()
