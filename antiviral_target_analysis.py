#!/usr/bin/env python3
"""
HMPV Antiviral Target Identification
=====================================

This script performs systematic gene and reaction knockout analysis to identify
potential antiviral targets for Human Metapneumovirus (HMPV). It identifies host 
genes and reactions that are essential for viral production.

Methodology:
------------
1. Load the integrated host-virus metabolic model
2. Perform single gene knockout analysis for all host genes
3. Perform single reaction knockout analysis for all reactions
4. Identify essential genes/reactions (lethal or significant impact on VBOF)
5. Analyze metabolic subsystem essentiality
6. Generate comprehensive reports

Output Files:
-------------
- gene_knockout_results.csv: Complete gene knockout analysis
- reaction_knockout_results.csv: Complete reaction knockout analysis
- top_gene_targets.csv: Top antiviral gene targets
- top_reaction_targets.csv: Top antiviral reaction targets
- subsystem_essentiality.csv: Subsystem essentiality analysis
- antiviral_targets_report.txt: Comprehensive report

Usage:
------
    python antiviral_target_analysis.py

Author: Syed Mushahid Hussain
"""

import logging
import sys
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import warnings

import cobra
from cobra import Model
import pandas as pd
import numpy as np

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONSTANTS
# ============================================================================
VBOF_REACTION_ID = "HMPV_VBOF"
ESSENTIALITY_THRESHOLD = 0.01      # < 1% of max flux = lethal
SIGNIFICANT_THRESHOLD = 0.5         # < 50% of max flux = significant
MODERATE_THRESHOLD = 0.9            # < 90% of max flux = moderate


# ============================================================================
# MODEL LOADING
# ============================================================================
def load_integrated_model(model_path: Path) -> Model:
    """
    Load the integrated host model with HMPV VBOF.
    
    Parameters:
    -----------
    model_path : Path
        Path to the integrated model SBML file.
    
    Returns:
    --------
    cobra.Model : Loaded model with VBOF as objective
    
    Raises:
    -------
    ValueError : If VBOF reaction not found
    FileNotFoundError : If model file not found
    """
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    logger.info(f"Loading model from: {model_path}")
    
    # Suppress COBRApy warning about missing objective (we set it below)
    cobra_logger = logging.getLogger('cobra.io.sbml')
    original_level = cobra_logger.level
    cobra_logger.setLevel(logging.CRITICAL)
    
    try:
        model = cobra.io.read_sbml_model(str(model_path))
    finally:
        cobra_logger.setLevel(original_level)
    
    if VBOF_REACTION_ID not in model.reactions:
        raise ValueError(f"VBOF reaction '{VBOF_REACTION_ID}' not found in model")
    
    model.objective = VBOF_REACTION_ID
    logger.info(f"Objective set to: {VBOF_REACTION_ID}")
    logger.info(f"Model loaded: {len(model.reactions)} reactions, {len(model.genes)} genes")
    
    return model


def get_baseline_flux(model: Model) -> float:
    """
    Calculate baseline VBOF flux (wild-type).
    
    Parameters:
    -----------
    model : Model
        The metabolic model
    
    Returns:
    --------
    float : Maximum VBOF flux (wild-type)
    """
    solution = model.optimize()
    if solution.status == 'optimal':
        return solution.objective_value
    return 0.0


# ============================================================================
# GENE KNOCKOUT ANALYSIS
# ============================================================================
def perform_gene_knockout_analysis(
    model: Model, 
    baseline_flux: float
) -> pd.DataFrame:
    """
    Perform single gene knockout analysis for all genes.
    
    Parameters:
    -----------
    model : Model
        The metabolic model
    baseline_flux : float
        Wild-type VBOF flux
    
    Returns:
    --------
    pd.DataFrame : Gene knockout results with columns:
        - gene_id, gene_name, knockout_flux, baseline_flux
        - flux_ratio, flux_reduction, impact, status
        - num_reactions, reactions
    """
    logger.info("Starting single gene knockout analysis...")
    logger.info(f"Total genes in model: {len(model.genes)}")
    
    knockout_data = []
    total_genes = len(model.genes)
    
    for i, gene in enumerate(model.genes):
        if (i + 1) % 500 == 0:
            logger.info(f"Progress: {i+1}/{total_genes} genes tested...")
        
        # Perform knockout using context manager
        with model:
            gene.knock_out()
            try:
                solution = model.optimize()
                flux = solution.objective_value if solution.status == 'optimal' else 0.0
                status = solution.status
            except Exception:
                flux = 0.0
                status = 'error'
        
        # Calculate metrics
        flux_ratio = flux / baseline_flux if baseline_flux > 0 else 0
        flux_reduction = 1 - flux_ratio
        
        # Classify impact
        if flux_ratio < ESSENTIALITY_THRESHOLD:
            impact = 'LETHAL'
        elif flux_ratio < SIGNIFICANT_THRESHOLD:
            impact = 'SIGNIFICANT'
        elif flux_ratio < MODERATE_THRESHOLD:
            impact = 'MODERATE'
        else:
            impact = 'MINIMAL'
        
        # Get gene info
        associated_rxns = [r.id for r in gene.reactions]
        
        knockout_data.append({
            'gene_id': gene.id,
            'gene_name': gene.name if gene.name else gene.id,
            'knockout_flux': flux,
            'baseline_flux': baseline_flux,
            'flux_ratio': flux_ratio,
            'flux_reduction': flux_reduction,
            'impact': impact,
            'status': status,
            'num_reactions': len(associated_rxns),
            'reactions': ';'.join(associated_rxns[:10])
        })
    
    df = pd.DataFrame(knockout_data)
    df = df.sort_values('flux_reduction', ascending=False)
    
    # Log summary
    logger.info("Gene knockout analysis complete!")
    logger.info(f"  Lethal knockouts: {len(df[df['impact'] == 'LETHAL'])}")
    logger.info(f"  Significant impact: {len(df[df['impact'] == 'SIGNIFICANT'])}")
    logger.info(f"  Moderate impact: {len(df[df['impact'] == 'MODERATE'])}")
    logger.info(f"  Minimal impact: {len(df[df['impact'] == 'MINIMAL'])}")
    
    return df


# ============================================================================
# REACTION KNOCKOUT ANALYSIS
# ============================================================================
def perform_reaction_knockout_analysis(
    model: Model, 
    baseline_flux: float
) -> pd.DataFrame:
    """
    Perform single reaction knockout analysis for all reactions.
    
    Parameters:
    -----------
    model : Model
        The metabolic model
    baseline_flux : float
        Wild-type VBOF flux
    
    Returns:
    --------
    pd.DataFrame : Reaction knockout results with columns:
        - reaction_id, reaction_name, subsystem
        - knockout_flux, baseline_flux, flux_ratio, flux_reduction
        - impact, status, num_genes, gene_ids, gene_names, equation
    """
    logger.info("Starting single reaction knockout analysis...")
    
    # Exclude exchange reactions and VBOF itself
    reactions_to_test = [r for r in model.reactions 
                        if not r.id.startswith('EX_') 
                        and r.id != VBOF_REACTION_ID]
    
    logger.info(f"Testing {len(reactions_to_test)} reactions...")
    
    knockout_data = []
    total_rxns = len(reactions_to_test)
    
    for i, rxn in enumerate(reactions_to_test):
        if (i + 1) % 500 == 0:
            logger.info(f"Progress: {i+1}/{total_rxns} reactions tested...")
        
        # Perform knockout using context manager
        with model:
            rxn.knock_out()
            try:
                solution = model.optimize()
                flux = solution.objective_value if solution.status == 'optimal' else 0.0
                status = solution.status
            except Exception:
                flux = 0.0
                status = 'error'
        
        # Calculate metrics
        flux_ratio = flux / baseline_flux if baseline_flux > 0 else 0
        flux_reduction = 1 - flux_ratio
        
        # Classify impact
        if flux_ratio < ESSENTIALITY_THRESHOLD:
            impact = 'LETHAL'
        elif flux_ratio < SIGNIFICANT_THRESHOLD:
            impact = 'SIGNIFICANT'
        elif flux_ratio < MODERATE_THRESHOLD:
            impact = 'MODERATE'
        else:
            impact = 'MINIMAL'
        
        # Get reaction info
        associated_genes = [g.id for g in rxn.genes]
        gene_names = [g.name if g.name else g.id for g in rxn.genes]
        
        knockout_data.append({
            'reaction_id': rxn.id,
            'reaction_name': rxn.name if rxn.name else rxn.id,
            'subsystem': rxn.subsystem if rxn.subsystem else 'Unknown',
            'knockout_flux': flux,
            'baseline_flux': baseline_flux,
            'flux_ratio': flux_ratio,
            'flux_reduction': flux_reduction,
            'impact': impact,
            'status': status,
            'num_genes': len(associated_genes),
            'gene_ids': ';'.join(associated_genes[:10]),
            'gene_names': ';'.join(gene_names[:10]),
            'equation': rxn.reaction
        })
    
    df = pd.DataFrame(knockout_data)
    df = df.sort_values('flux_reduction', ascending=False)
    
    # Log summary
    logger.info("Reaction knockout analysis complete!")
    logger.info(f"  Lethal knockouts: {len(df[df['impact'] == 'LETHAL'])}")
    logger.info(f"  Significant impact: {len(df[df['impact'] == 'SIGNIFICANT'])}")
    logger.info(f"  Moderate impact: {len(df[df['impact'] == 'MODERATE'])}")
    logger.info(f"  Minimal impact: {len(df[df['impact'] == 'MINIMAL'])}")
    
    return df


# ============================================================================
# SUBSYSTEM ANALYSIS
# ============================================================================
def analyze_subsystem_essentiality(reaction_results: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze which metabolic subsystems are most essential for viral production.
    
    Parameters:
    -----------
    reaction_results : pd.DataFrame
        Reaction knockout results
    
    Returns:
    --------
    pd.DataFrame : Subsystem essentiality summary
    """
    logger.info("Analyzing subsystem essentiality...")
    
    subsystem_data = []
    
    for subsystem in reaction_results['subsystem'].unique():
        subsystem_rxns = reaction_results[reaction_results['subsystem'] == subsystem]
        
        lethal_count = len(subsystem_rxns[subsystem_rxns['impact'] == 'LETHAL'])
        significant_count = len(subsystem_rxns[subsystem_rxns['impact'] == 'SIGNIFICANT'])
        total_count = len(subsystem_rxns)
        
        avg_reduction = subsystem_rxns['flux_reduction'].mean()
        max_reduction = subsystem_rxns['flux_reduction'].max()
        
        # Get sample lethal reactions
        lethal_rxns = subsystem_rxns[subsystem_rxns['impact'] == 'LETHAL']['reaction_name'].head(3).tolist()
        
        essentiality_score = (lethal_count * 2 + significant_count) / total_count if total_count > 0 else 0
        
        subsystem_data.append({
            'subsystem': subsystem,
            'total_reactions': total_count,
            'lethal_knockouts': lethal_count,
            'significant_knockouts': significant_count,
            'avg_flux_reduction': avg_reduction,
            'max_flux_reduction': max_reduction,
            'essentiality_score': essentiality_score,
            'sample_lethal_reactions': ';'.join(lethal_rxns)
        })
    
    df = pd.DataFrame(subsystem_data)
    df = df.sort_values('essentiality_score', ascending=False)
    
    return df


# ============================================================================
# TOP TARGETS EXTRACTION
# ============================================================================
def extract_top_targets(
    gene_results: pd.DataFrame,
    reaction_results: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract top antiviral targets from knockout results.
    
    Parameters:
    -----------
    gene_results : pd.DataFrame
        Gene knockout results
    reaction_results : pd.DataFrame
        Reaction knockout results
    
    Returns:
    --------
    Tuple[pd.DataFrame, pd.DataFrame] : Top gene and reaction targets
    """
    logger.info("Extracting top antiviral targets...")
    
    # Get lethal and significant gene targets
    lethal_genes = gene_results[gene_results['impact'] == 'LETHAL'].copy()
    significant_genes = gene_results[gene_results['impact'] == 'SIGNIFICANT'].copy()
    top_genes = pd.concat([lethal_genes, significant_genes]).sort_values(
        'flux_reduction', ascending=False
    )
    
    # Get lethal and significant reaction targets
    lethal_rxns = reaction_results[reaction_results['impact'] == 'LETHAL'].copy()
    significant_rxns = reaction_results[reaction_results['impact'] == 'SIGNIFICANT'].copy()
    top_rxns = pd.concat([lethal_rxns, significant_rxns]).sort_values(
        'flux_reduction', ascending=False
    )
    
    logger.info(f"Top gene targets: {len(top_genes)}")
    logger.info(f"Top reaction targets: {len(top_rxns)}")
    
    return top_genes, top_rxns


# ============================================================================
# REPORT GENERATION
# ============================================================================
def generate_report(
    model: Model,
    baseline_flux: float,
    gene_results: pd.DataFrame,
    reaction_results: pd.DataFrame,
    top_genes: pd.DataFrame,
    top_rxns: pd.DataFrame,
    subsystem_analysis: pd.DataFrame,
    output_path: Path
) -> str:
    """
    Generate comprehensive antiviral target report.
    
    Parameters:
    -----------
    Various analysis results and output path.
    
    Returns:
    --------
    str : Report text
    """
    report = f"""
================================================================================
HMPV ANTIVIRAL TARGET IDENTIFICATION REPORT
================================================================================

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MODEL INFORMATION:
------------------
Model ID: {model.id}
Total Reactions: {len(model.reactions)}
Total Genes: {len(model.genes)}
VBOF Reaction: {VBOF_REACTION_ID}
Baseline VBOF Flux: {baseline_flux:.6f}

================================================================================
SUMMARY
================================================================================

GENE KNOCKOUTS:
  Total tested: {len(gene_results)}
  Lethal (>99% reduction): {len(gene_results[gene_results['impact'] == 'LETHAL'])}
  Significant (50-99%): {len(gene_results[gene_results['impact'] == 'SIGNIFICANT'])}
  Moderate (10-50%): {len(gene_results[gene_results['impact'] == 'MODERATE'])}
  Minimal (<10%): {len(gene_results[gene_results['impact'] == 'MINIMAL'])}

REACTION KNOCKOUTS:
  Total tested: {len(reaction_results)}
  Lethal (>99% reduction): {len(reaction_results[reaction_results['impact'] == 'LETHAL'])}
  Significant (50-99%): {len(reaction_results[reaction_results['impact'] == 'SIGNIFICANT'])}
  Moderate (10-50%): {len(reaction_results[reaction_results['impact'] == 'MODERATE'])}
  Minimal (<10%): {len(reaction_results[reaction_results['impact'] == 'MINIMAL'])}

================================================================================
ESSENTIAL GENE TARGETS (Lethal for HMPV Production)
================================================================================
"""
    
    for idx, (i, row) in enumerate(top_genes[top_genes['impact'] == 'LETHAL'].head(20).iterrows()):
        report += f"""
{idx + 1}. {row['gene_id']} - {row['gene_name']}
   Flux Reduction: {row['flux_reduction']*100:.1f}%
   Associated Reactions ({row['num_reactions']}): {row['reactions']}
"""

    report += """
================================================================================
ESSENTIAL REACTION TARGETS (Lethal for HMPV Production)
================================================================================
"""
    
    for idx, (i, row) in enumerate(top_rxns[top_rxns['impact'] == 'LETHAL'].head(20).iterrows()):
        equation_short = row['equation'][:80] + '...' if len(str(row['equation'])) > 80 else row['equation']
        report += f"""
{idx + 1}. {row['reaction_id']} - {row['reaction_name']}
   Subsystem: {row['subsystem']}
   Flux Reduction: {row['flux_reduction']*100:.1f}%
   Genes: {row['gene_names']}
   Equation: {equation_short}
"""

    # Add significant targets if any exist
    significant_genes = top_genes[top_genes['impact'] == 'SIGNIFICANT']
    if len(significant_genes) > 0:
        report += """
================================================================================
SIGNIFICANT IMPACT GENE TARGETS (50-99% reduction)
================================================================================
"""
        for idx, (i, row) in enumerate(significant_genes.head(15).iterrows()):
            report += f"""
{idx + 1}. {row['gene_id']} - {row['gene_name']}
   Flux Reduction: {row['flux_reduction']*100:.1f}%
   Associated Reactions: {row['reactions']}
"""

    significant_rxns = top_rxns[top_rxns['impact'] == 'SIGNIFICANT']
    if len(significant_rxns) > 0:
        report += """
================================================================================
SIGNIFICANT IMPACT REACTION TARGETS (50-99% reduction)
================================================================================
"""
        for idx, (i, row) in enumerate(significant_rxns.head(15).iterrows()):
            report += f"""
{idx + 1}. {row['reaction_id']} - {row['reaction_name']}
   Flux Reduction: {row['flux_reduction']*100:.1f}%
   Genes: {row['gene_names']}
"""

    # Subsystem analysis
    essential_subsystems = subsystem_analysis[subsystem_analysis['lethal_knockouts'] > 0]
    if len(essential_subsystems) > 0:
        report += """
================================================================================
SUBSYSTEM ESSENTIALITY ANALYSIS
================================================================================
"""
        for idx, (i, row) in enumerate(essential_subsystems.head(15).iterrows()):
            report += f"""
{idx + 1}. {row['subsystem']}
   Total Reactions: {row['total_reactions']}
   Lethal Knockouts: {row['lethal_knockouts']}
   Essentiality Score: {row['essentiality_score']:.3f}
   Sample Lethal Reactions: {row['sample_lethal_reactions']}
"""

    report += f"""
================================================================================
INTERPRETATION AND RECOMMENDATIONS
================================================================================

The analysis identified {len(top_genes[top_genes['impact'] == 'LETHAL'])} essential host genes and 
{len(top_rxns[top_rxns['impact'] == 'LETHAL'])} essential reactions for HMPV virion production.

POTENTIAL DRUG TARGETS:
-----------------------
Genes and reactions with lethal or significant impact on viral production
represent potential antiviral drug targets. Priority targets are those that:
1. Show >50% reduction in viral production when knocked out
2. Have known inhibitors or are druggable
3. Are not essential for host cell survival


================================================================================
FILES GENERATED
================================================================================
1. gene_knockout_results.csv - Complete gene knockout analysis
2. reaction_knockout_results.csv - Complete reaction knockout analysis
3. top_gene_targets.csv - Top antiviral gene targets
4. top_reaction_targets.csv - Top antiviral reaction targets
5. subsystem_essentiality.csv - Subsystem essentiality analysis
6. antiviral_targets_report.txt - This report

================================================================================
"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    return report


# ============================================================================
# MAIN PIPELINE
# ============================================================================
def main():
    """
    Main pipeline for antiviral target identification.
    """
    logger.info("=" * 70)
    logger.info("HMPV ANTIVIRAL TARGET IDENTIFICATION")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Define paths
    model_path = Path("output/iHsaEC21_CLEAN_with_HMPV_VBOF.xml")
    output_dir = Path("output/antiviral_analysis")
    output_dir.mkdir(exist_ok=True)
    
    # =========================================================================
    # STEP 1: Load model
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 1: LOADING MODEL")
    logger.info("=" * 70)
    
    model = load_integrated_model(model_path)
    
    # =========================================================================
    # STEP 2: Calculate baseline flux
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: CALCULATING BASELINE FLUX")
    logger.info("=" * 70)
    
    baseline_flux = get_baseline_flux(model)
    logger.info(f"Baseline VBOF flux (wild-type): {baseline_flux:.6f}")
    
    if baseline_flux <= 0:
        logger.error("Baseline flux is zero or negative. Cannot proceed.")
        sys.exit(1)
    
    # =========================================================================
    # STEP 3: Gene knockout analysis
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: GENE KNOCKOUT ANALYSIS")
    logger.info("=" * 70)
    
    gene_results = perform_gene_knockout_analysis(model, baseline_flux)
    gene_results.to_csv(output_dir / "gene_knockout_results.csv", index=False)
    logger.info(f"Saved: {output_dir / 'gene_knockout_results.csv'}")
    
    # =========================================================================
    # STEP 4: Reaction knockout analysis
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 4: REACTION KNOCKOUT ANALYSIS")
    logger.info("=" * 70)
    
    reaction_results = perform_reaction_knockout_analysis(model, baseline_flux)
    reaction_results.to_csv(output_dir / "reaction_knockout_results.csv", index=False)
    logger.info(f"Saved: {output_dir / 'reaction_knockout_results.csv'}")
    
    # =========================================================================
    # STEP 5: Extract top targets
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 5: EXTRACTING TOP TARGETS")
    logger.info("=" * 70)
    
    top_genes, top_rxns = extract_top_targets(gene_results, reaction_results)
    top_genes.to_csv(output_dir / "top_gene_targets.csv", index=False)
    top_rxns.to_csv(output_dir / "top_reaction_targets.csv", index=False)
    logger.info(f"Saved: {output_dir / 'top_gene_targets.csv'}")
    logger.info(f"Saved: {output_dir / 'top_reaction_targets.csv'}")
    
    # =========================================================================
    # STEP 6: Subsystem analysis
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 6: SUBSYSTEM ANALYSIS")
    logger.info("=" * 70)
    
    subsystem_analysis = analyze_subsystem_essentiality(reaction_results)
    subsystem_analysis.to_csv(output_dir / "subsystem_essentiality.csv", index=False)
    logger.info(f"Saved: {output_dir / 'subsystem_essentiality.csv'}")
    
    # =========================================================================
    # STEP 7: Generate report
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 7: GENERATING REPORT")
    logger.info("=" * 70)
    
    report = generate_report(
        model, baseline_flux, gene_results, reaction_results,
        top_genes, top_rxns, subsystem_analysis,
        output_dir / "antiviral_targets_report.txt"
    )
    logger.info(f"Saved: {output_dir / 'antiviral_targets_report.txt'}")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    lethal_genes = gene_results[gene_results['impact'] == 'LETHAL']
    lethal_rxns = reaction_results[reaction_results['impact'] == 'LETHAL']
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nBaseline VBOF Flux: {baseline_flux:.6f}")
    
    print(f"\n--- LETHAL GENE TARGETS ({len(lethal_genes)}) ---")
    for _, row in lethal_genes.head(10).iterrows():
        print(f"  {row['gene_id']} ({row['gene_name']})")
    
    print(f"\n--- LETHAL REACTION TARGETS ({len(lethal_rxns)}) ---")
    for _, row in lethal_rxns.head(10).iterrows():
        print(f"  {row['reaction_id']} ({row['reaction_name']})")
    
    print(f"\nOutput saved to: {output_dir}")
    print("=" * 70)
    
    logger.info("\n" + "=" * 70)
    logger.info("ANTIVIRAL TARGET ANALYSIS COMPLETE")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
