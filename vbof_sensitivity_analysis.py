#!/usr/bin/env python3
"""
HMPV VBOF Sensitivity Analysis
================================

This script performs sensitivity analysis by testing different VBOF parameter
combinations and comparing antiviral targets across scenarios.

Key Parameters Varied:
- Virion diameter (nm): Affects lipid envelope requirements
- F protein copy number: Affects glycan requirements
- G protein copy number: Affects glycan requirements

For each parameter combination:
1. Build VBOF with specified parameters
2. Integrate into host model
3. Run antiviral target analysis
4. Compare results across scenarios

Output:
-------
- output/sensitivity_analysis/: Directory with results for each scenario
- output/sensitivity_analysis/comparison_report.txt: Cross-scenario comparison
- output/sensitivity_analysis/robust_targets.csv: Targets appearing in multiple scenarios
- output/sensitivity_analysis/scenario_summary.csv: Summary of all scenarios

Usage:
------
    python vbof_sensitivity_analysis.py

Author: Syed Mushahid Hussain
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from itertools import product
import warnings

import pandas as pd
import numpy as np
from cobra import Model

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.genome_analyzer import (
    load_genome,
    parse_gff_annotations,
    calculate_genome_stoichiometry,
    get_genome_summary
)
from src.protein_analyzer import (
    load_proteins,
    calculate_protein_stoichiometry,
    get_protein_summary
)
from src.vbof_builder import (
    build_vbof,
    export_vbof_to_dict
)
from src.model_integration import (
    load_host_model,
    integrate_vbof,
    save_integrated_model
)
from src.config import (
    GENOMIC_DIR,
    PROTEIN_DIR,
    MODEL_DIR,
    OUTPUT_DIR,
    DEFAULT_GENOME_FILE,
    DEFAULT_GFF_FILE,
    DEFAULT_PROTEIN_FILE,
    DEFAULT_HOST_MODEL,
    HMPV_COPY_NUMBERS,
    VBOF_REACTION_ID
)
from antiviral_target_analysis import (
    load_integrated_model,
    get_baseline_flux,
    perform_gene_knockout_analysis,
    perform_reaction_knockout_analysis,
    extract_top_targets
)

# Suppress warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# ============================================================================
# PARAMETER RANGES FOR SENSITIVITY ANALYSIS
# ============================================================================

# Virion diameter range (nm)
# HMPV is pleomorphic: 150-300 nm typical, up to 600 nm for filamentous
VIRION_DIAMETER_RANGE = [150.0, 200.0, 250.0, 300.0]

# F protein copy number range
# Based on cryo-EM studies of paramyxoviruses
F_COPY_NUMBER_RANGE = [250, 350, 450]

# G protein copy number range
# G protein is variable and can be absent in some strains
G_COPY_NUMBER_RANGE = [150, 250, 350]

# ============================================================================
# SCENARIO MANAGEMENT
# ============================================================================

class Scenario:
    """Container for a single VBOF parameter scenario."""
    
    def __init__(
        self,
        scenario_id: str,
        virion_diameter_nm: float,
        f_copy_number: int,
        g_copy_number: int
    ):
        self.scenario_id = scenario_id
        self.virion_diameter_nm = virion_diameter_nm
        self.f_copy_number = f_copy_number
        self.g_copy_number = g_copy_number
        self.vbof = None
        self.model = None
        self.baseline_flux = None
        self.gene_results = None
        self.reaction_results = None
        self.top_genes = None
        self.top_rxns = None
    
    def __repr__(self):
        return (f"Scenario({self.scenario_id}: "
                f"d={self.virion_diameter_nm}nm, "
                f"F={self.f_copy_number}, G={self.g_copy_number})")


def generate_scenarios() -> List[Scenario]:
    """
    Generate all parameter combinations for sensitivity analysis.
    
    Returns:
    --------
    List[Scenario] : List of scenario objects
    """
    scenarios = []
    scenario_num = 1
    
    for diameter, f_copy, g_copy in product(
        VIRION_DIAMETER_RANGE,
        F_COPY_NUMBER_RANGE,
        G_COPY_NUMBER_RANGE
    ):
        scenario_id = f"scenario_{scenario_num:03d}"
        scenario = Scenario(
            scenario_id=scenario_id,
            virion_diameter_nm=diameter,
            f_copy_number=f_copy,
            g_copy_number=g_copy
        )
        scenarios.append(scenario)
        scenario_num += 1
    
    logger.info(f"Generated {len(scenarios)} scenarios for sensitivity analysis")
    return scenarios


# ============================================================================
# VBOF BUILDING AND INTEGRATION
# ============================================================================

def build_vbof_for_scenario(
    scenario: Scenario,
    genome_stoichiometry: Dict[str, float],
    protein_stoichiometry: Dict[str, float],
    genome_length: int,
    total_amino_acids: int,
    num_proteins: int,
    base_copy_numbers: Dict[str, int]
) -> None:
    """
    Build VBOF for a specific scenario.
    
    Parameters:
    -----------
    scenario : Scenario
        Scenario object to populate
    genome_stoichiometry : dict
        Genome nucleotide stoichiometry (constant across scenarios)
    protein_stoichiometry : dict
        Base protein stoichiometry (will be adjusted for F/G copy numbers)
    genome_length : int
        Genome length
    total_amino_acids : int
        Total amino acids (will be recalculated)
    num_proteins : int
        Number of proteins
    base_copy_numbers : dict
        Base copy numbers for all proteins
    """
    logger.info(f"Building VBOF for {scenario.scenario_id}...")
    logger.info(f"  Parameters: d={scenario.virion_diameter_nm}nm, "
                f"F={scenario.f_copy_number}, G={scenario.g_copy_number}")
    
    # Adjust copy numbers for F and G
    adjusted_copy_numbers = base_copy_numbers.copy()
    adjusted_copy_numbers['F'] = scenario.f_copy_number
    adjusted_copy_numbers['G'] = scenario.g_copy_number
    
    # Recalculate protein stoichiometry with adjusted copy numbers
    # Load proteins to recalculate properly
    protein_path = PROTEIN_DIR / DEFAULT_PROTEIN_FILE
    proteins = load_proteins(str(protein_path))
    
    # Recalculate protein stoichiometry with new copy numbers
    adjusted_protein_stoichiometry = calculate_protein_stoichiometry(
        proteins,
        copy_numbers=adjusted_copy_numbers
    )
    
    # Recalculate total amino acids
    adjusted_total_aa = sum(
        proteins[gene].length * adjusted_copy_numbers[gene]
        for gene in proteins.keys()
        if gene in adjusted_copy_numbers
    )
    
    # Build VBOF
    scenario.vbof = build_vbof(
        genome_stoichiometry=genome_stoichiometry,
        protein_stoichiometry=adjusted_protein_stoichiometry,
        genome_length=genome_length,
        total_amino_acids=adjusted_total_aa,
        num_proteins=num_proteins,
        virion_diameter_nm=scenario.virion_diameter_nm,
        include_lipids=True,
        include_energy=True,
        include_glycans=True,
        f_copy_number=scenario.f_copy_number,
        g_copy_number=scenario.g_copy_number
    )
    
    logger.info(f"  VBOF built: {scenario.vbof.metadata['total_metabolites']} metabolites")


def normalize_vbof_stoichiometry(vbof_stoichiometry: Dict[str, float]) -> Dict[str, float]:
    """
    Normalize VBOF stoichiometry using simple normalization.
    
    Parameters:
    -----------
    vbof_stoichiometry : dict
        Raw VBOF stoichiometry
    
    Returns:
    --------
    dict : Normalized stoichiometry
    """
    # Calculate total consumed for normalization
    total_consumed = sum(abs(c) for c in vbof_stoichiometry.values() if c < 0)
    
    if total_consumed == 0:
        logger.warning("  Total consumed is zero, cannot normalize")
        return vbof_stoichiometry
    
    # Normalize by dividing by total consumed
    normalized = {
        met_id: coef / total_consumed
        for met_id, coef in vbof_stoichiometry.items()
    }
    
    return normalized


def integrate_vbof_for_scenario(
    scenario: Scenario,
    model_path: Path
) -> None:
    """
    Integrate VBOF into host model for a scenario.
    
    Parameters:
    -----------
    scenario : Scenario
        Scenario with built VBOF
    model_path : Path
        Path to host model file (will be loaded fresh for each scenario)
    """
    logger.info(f"Integrating VBOF for {scenario.scenario_id}...")
    
    # Load fresh model copy to avoid deepcopy issues with GPR objects
    model_copy = load_host_model(model_path)
    
    # Normalize VBOF before integration
    raw_stoichiometry = scenario.vbof.combined_stoichiometry
    vbof_stoichiometry = normalize_vbof_stoichiometry(raw_stoichiometry)
    
    logger.info(f"  Normalized VBOF (total consumed: {sum(abs(c) for c in vbof_stoichiometry.values() if c < 0):.6f})")
    
    # Integrate VBOF
    model_copy, vbof_reaction, unmapped = integrate_vbof(
        model=model_copy,
        vbof_stoichiometry=vbof_stoichiometry,
        reaction_id=VBOF_REACTION_ID,
        skip_unmapped=True
    )
    
    # Set objective
    model_copy.objective = VBOF_REACTION_ID
    
    scenario.model = model_copy
    
    if unmapped:
        logger.warning(f"  {len(unmapped)} metabolites unmapped: {unmapped[:5]}...")


# ============================================================================
# ANTIVIRAL ANALYSIS
# ============================================================================

def run_antiviral_analysis_for_scenario(
    scenario: Scenario,
    output_dir: Path
) -> None:
    """
    Run antiviral target analysis for a scenario.
    
    Parameters:
    -----------
    scenario : Scenario
        Scenario with integrated model
    output_dir : Path
        Output directory for this scenario
    """
    logger.info(f"Running antiviral analysis for {scenario.scenario_id}...")
    
    # Calculate baseline flux
    try:
        scenario.baseline_flux = get_baseline_flux(scenario.model)
        logger.info(f"  Baseline flux: {scenario.baseline_flux:.6f}")
        
        if scenario.baseline_flux <= 0:
            logger.warning(f"  Baseline flux is zero or negative! Status may be non-optimal.")
            # Try to get solution status
            solution = scenario.model.optimize()
            logger.warning(f"  Optimization status: {solution.status}")
            if solution.status != 'optimal':
                logger.error(f"  Model optimization failed with status: {solution.status}. Skipping analysis.")
                return
            else:
                logger.warning(f"  Model is optimal but flux is zero. This may indicate VBOF is infeasible.")
                return
    except Exception as e:
        logger.error(f"  Error calculating baseline flux: {e}")
        return
    
    # Run gene knockout analysis
    logger.info("  Running gene knockout analysis...")
    scenario.gene_results = perform_gene_knockout_analysis(
        scenario.model,
        scenario.baseline_flux
    )
    
    # Run reaction knockout analysis
    logger.info("  Running reaction knockout analysis...")
    scenario.reaction_results = perform_reaction_knockout_analysis(
        scenario.model,
        scenario.baseline_flux
    )
    
    # Extract top targets
    scenario.top_genes, scenario.top_rxns = extract_top_targets(
        scenario.gene_results,
        scenario.reaction_results
    )
    
    # Save results
    output_dir.mkdir(exist_ok=True, parents=True)
    scenario.gene_results.to_csv(output_dir / "gene_knockout_results.csv", index=False)
    scenario.reaction_results.to_csv(output_dir / "reaction_knockout_results.csv", index=False)
    scenario.top_genes.to_csv(output_dir / "top_gene_targets.csv", index=False)
    scenario.top_rxns.to_csv(output_dir / "top_reaction_targets.csv", index=False)
    
    logger.info(f"  Results saved to: {output_dir}")


# ============================================================================
# COMPARISON AND AGGREGATION
# ============================================================================

def compare_scenarios(scenarios: List[Scenario]) -> Dict:
    """
    Compare results across all scenarios to identify robust targets.
    
    Parameters:
    -----------
    scenarios : List[Scenario]
        List of completed scenarios
    
    Returns:
    --------
    dict : Comparison results including robust targets
    """
    logger.info("Comparing results across scenarios...")
    
    # Collect all lethal targets
    all_lethal_genes = set()
    all_lethal_rxns = set()
    all_significant_genes = set()
    all_significant_rxns = set()
    
    scenario_gene_counts = {}
    scenario_rxn_counts = {}
    
    for scenario in scenarios:
        if scenario.top_genes is None or scenario.top_rxns is None:
            continue
        
        # Collect lethal genes
        lethal_genes = set(scenario.top_genes[scenario.top_genes['impact'] == 'LETHAL']['gene_id'])
        all_lethal_genes.update(lethal_genes)
        for gene in lethal_genes:
            scenario_gene_counts[gene] = scenario_gene_counts.get(gene, 0) + 1
        
        # Collect lethal reactions
        lethal_rxns = set(scenario.top_rxns[scenario.top_rxns['impact'] == 'LETHAL']['reaction_id'])
        all_lethal_rxns.update(lethal_rxns)
        for rxn in lethal_rxns:
            scenario_rxn_counts[rxn] = scenario_rxn_counts.get(rxn, 0) + 1
        
        # Collect significant targets
        sig_genes = set(scenario.top_genes[scenario.top_genes['impact'] == 'SIGNIFICANT']['gene_id'])
        all_significant_genes.update(sig_genes)
        
        sig_rxns = set(scenario.top_rxns[scenario.top_rxns['impact'] == 'SIGNIFICANT']['reaction_id'])
        all_significant_rxns.update(sig_rxns)
    
    num_scenarios = len([s for s in scenarios if s.top_genes is not None])
    
    # Calculate robustness scores (fraction of scenarios where target is lethal)
    robust_genes = []
    for gene in all_lethal_genes:
        count = scenario_gene_counts.get(gene, 0)
        robustness = count / num_scenarios if num_scenarios > 0 else 0
        robust_genes.append({
            'gene_id': gene,
            'scenarios_lethal': count,
            'robustness_score': robustness,
            'appears_in_all': count == num_scenarios
        })
    
    robust_rxns = []
    for rxn in all_lethal_rxns:
        count = scenario_rxn_counts.get(rxn, 0)
        robustness = count / num_scenarios if num_scenarios > 0 else 0
        robust_rxns.append({
            'reaction_id': rxn,
            'scenarios_lethal': count,
            'robustness_score': robustness,
            'appears_in_all': count == num_scenarios
        })
    
    # Create DataFrames with proper columns even if empty
    if robust_genes:
        robust_genes_df = pd.DataFrame(robust_genes)
        robust_genes_df = robust_genes_df.sort_values('robustness_score', ascending=False)
    else:
        robust_genes_df = pd.DataFrame(columns=['gene_id', 'scenarios_lethal', 'robustness_score', 'appears_in_all'])
    
    if robust_rxns:
        robust_rxns_df = pd.DataFrame(robust_rxns)
        robust_rxns_df = robust_rxns_df.sort_values('robustness_score', ascending=False)
    else:
        robust_rxns_df = pd.DataFrame(columns=['reaction_id', 'scenarios_lethal', 'robustness_score', 'appears_in_all'])
    
    return {
        'robust_genes': robust_genes_df,
        'robust_rxns': robust_rxns_df,
        'total_scenarios': num_scenarios,
        'total_lethal_genes': len(all_lethal_genes),
        'total_lethal_rxns': len(all_lethal_rxns),
        'universal_lethal_genes': len([g for g in robust_genes if g['appears_in_all']]),
        'universal_lethal_rxns': len([r for r in robust_rxns if r['appears_in_all']])
    }


def generate_comparison_report(
    scenarios: List[Scenario],
    comparison_results: Dict,
    output_path: Path
) -> None:
    """
    Generate comprehensive comparison report.
    
    Parameters:
    -----------
    scenarios : List[Scenario]
        All scenarios
    comparison_results : dict
        Comparison results from compare_scenarios
    output_path : Path
        Path to save report
    """
    report = f"""
================================================================================
HMPV VBOF SENSITIVITY ANALYSIS REPORT
================================================================================

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

================================================================================
PARAMETER RANGES TESTED
================================================================================

Virion Diameter (nm): {VIRION_DIAMETER_RANGE}
F Protein Copy Number: {F_COPY_NUMBER_RANGE}
G Protein Copy Number: {G_COPY_NUMBER_RANGE}

Total Scenarios: {len(scenarios)}
Completed Scenarios: {comparison_results['total_scenarios']}

================================================================================
SCENARIO SUMMARY
================================================================================
"""
    
    for scenario in scenarios:
        if scenario.baseline_flux is not None:
            report += f"""
{scenario.scenario_id}:
  Parameters: d={scenario.virion_diameter_nm}nm, F={scenario.f_copy_number}, G={scenario.g_copy_number}
  Baseline VBOF Flux: {scenario.baseline_flux:.6f}
  Lethal Genes: {len(scenario.top_genes[scenario.top_genes['impact'] == 'LETHAL']) if scenario.top_genes is not None else 'N/A'}
  Lethal Reactions: {len(scenario.top_rxns[scenario.top_rxns['impact'] == 'LETHAL']) if scenario.top_rxns is not None else 'N/A'}
"""
    
    report += f"""
================================================================================
ROBUST TARGETS ANALYSIS
================================================================================

Universal Lethal Targets (appear in ALL scenarios):
  Genes: {comparison_results['universal_lethal_genes']}
  Reactions: {comparison_results['universal_lethal_rxns']}

Total Unique Lethal Targets (across all scenarios):
  Genes: {comparison_results['total_lethal_genes']}
  Reactions: {comparison_results['total_lethal_rxns']}

================================================================================
TOP ROBUST GENE TARGETS
================================================================================

These genes are lethal in multiple scenarios, indicating they are robust
antiviral targets regardless of VBOF parameter uncertainty.

"""
    
    if len(comparison_results['robust_genes']) > 0:
        top_robust_genes = comparison_results['robust_genes'].head(30)
        for idx, row in top_robust_genes.iterrows():
            report += f"{row['gene_id']}: Lethal in {row['scenarios_lethal']}/{comparison_results['total_scenarios']} scenarios "
            report += f"(robustness: {row['robustness_score']*100:.1f}%)\n"
            if row['appears_in_all']:
                report += "  *** UNIVERSAL TARGET (appears in all scenarios) ***\n"
    else:
        report += "No robust gene targets found.\n"
    
    report += f"""
================================================================================
TOP ROBUST REACTION TARGETS
================================================================================

These reactions are lethal in multiple scenarios, indicating they are robust
antiviral targets regardless of VBOF parameter uncertainty.

"""
    
    if len(comparison_results['robust_rxns']) > 0:
        top_robust_rxns = comparison_results['robust_rxns'].head(30)
        for idx, row in top_robust_rxns.iterrows():
            report += f"{row['reaction_id']}: Lethal in {row['scenarios_lethal']}/{comparison_results['total_scenarios']} scenarios "
            report += f"(robustness: {row['robustness_score']*100:.1f}%)\n"
            if row['appears_in_all']:
                report += "  *** UNIVERSAL TARGET (appears in all scenarios) ***\n"
    else:
        report += "No robust reaction targets found.\n"
    
    report += f"""
================================================================================
INTERPRETATION
================================================================================

This sensitivity analysis tests the robustness of antiviral targets across
different VBOF parameter combinations. Targets that appear as lethal in
multiple scenarios are more reliable candidates for drug development, as they
are less sensitive to uncertainty in:

1. Virion size (pleomorphic nature of HMPV)
2. Protein copy numbers (experimental estimates)

RECOMMENDATIONS:
---------------
1. Prioritize universal targets (appear in all scenarios) for experimental validation
2. Consider robustness score when selecting targets for drug development
3. Targets with >70% robustness are likely reliable across parameter uncertainty
4. Further experimental data on virion size and copy numbers can refine these results

================================================================================
"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    logger.info(f"Comparison report saved to: {output_path}")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """
    Main sensitivity analysis pipeline.
    """
    logger.info("=" * 70)
    logger.info("HMPV VBOF SENSITIVITY ANALYSIS")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    output_dir = Path(OUTPUT_DIR) / "sensitivity_analysis"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    try:
        # =====================================================================
        # STEP 1: Load genome and protein data (constant across scenarios)
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 1: LOADING GENOME AND PROTEIN DATA")
        logger.info("=" * 70)
        
        genome_path = GENOMIC_DIR / DEFAULT_GENOME_FILE
        gff_path = GENOMIC_DIR / DEFAULT_GFF_FILE
        protein_path = PROTEIN_DIR / DEFAULT_PROTEIN_FILE
        
        # Load genome
        genome = load_genome(str(genome_path))
        annotations = parse_gff_annotations(str(gff_path))
        genome_stoichiometry = calculate_genome_stoichiometry(genome, copies_per_virion=1)
        genome_summary = get_genome_summary(genome)
        
        # Load proteins
        proteins = load_proteins(str(protein_path))
        protein_summary = get_protein_summary(proteins)
        
        # Calculate base protein stoichiometry
        base_protein_stoichiometry = calculate_protein_stoichiometry(
            proteins,
            copy_numbers=HMPV_COPY_NUMBERS
        )
        
        # Calculate base total amino acids
        base_total_aa = sum(
            proteins[gene].length * HMPV_COPY_NUMBERS[gene]
            for gene in proteins.keys()
            if gene in HMPV_COPY_NUMBERS
        )
        
        logger.info(f"Genome length: {genome.length} nt")
        logger.info(f"Total proteins: {len(proteins)}")
        logger.info(f"Base total amino acids: {base_total_aa}")
        
        # =====================================================================
        # STEP 2: Get host model path (will load fresh for each scenario)
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 2: PREPARING HOST MODEL")
        logger.info("=" * 70)
        
        model_path = MODEL_DIR / DEFAULT_HOST_MODEL
        # Load once to verify it exists and get info
        test_model = load_host_model(model_path)
        logger.info(f"Host model path: {model_path}")
        logger.info(f"Host model: {len(test_model.reactions)} reactions, "
                   f"{len(test_model.genes)} genes")
        # Don't keep the model in memory - will load fresh for each scenario
        
        # =====================================================================
        # STEP 3: Generate scenarios
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 3: GENERATING SCENARIOS")
        logger.info("=" * 70)
        
        scenarios = generate_scenarios()
        logger.info(f"Total scenarios to test: {len(scenarios)}")
        
        # =====================================================================
        # STEP 4: Process each scenario
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 4: PROCESSING SCENARIOS")
        logger.info("=" * 70)
        
        scenario_summaries = []
        
        for i, scenario in enumerate(scenarios, 1):
            logger.info(f"\n--- Processing {scenario.scenario_id} ({i}/{len(scenarios)}) ---")
            
            try:
                # Build VBOF
                build_vbof_for_scenario(
                    scenario=scenario,
                    genome_stoichiometry=genome_stoichiometry,
                    protein_stoichiometry=base_protein_stoichiometry,
                    genome_length=genome.length,
                    total_amino_acids=base_total_aa,
                    num_proteins=len(proteins),
                    base_copy_numbers=HMPV_COPY_NUMBERS
                )
                
                if scenario.vbof is None:
                    logger.error(f"  VBOF building failed for {scenario.scenario_id}")
                    continue
                
                # Integrate VBOF (loads fresh model copy)
                integrate_vbof_for_scenario(scenario, model_path)
                
                if scenario.model is None:
                    logger.error(f"  Model integration failed for {scenario.scenario_id}")
                    continue
                
                # Run antiviral analysis
                scenario_output_dir = output_dir / scenario.scenario_id
                run_antiviral_analysis_for_scenario(scenario, scenario_output_dir)
                
                # Save scenario summary (only if analysis completed successfully)
                if scenario.baseline_flux is not None and scenario.baseline_flux > 0:
                    if scenario.top_genes is not None and scenario.top_rxns is not None:
                        scenario_summaries.append({
                            'scenario_id': scenario.scenario_id,
                            'virion_diameter_nm': scenario.virion_diameter_nm,
                            'f_copy_number': scenario.f_copy_number,
                            'g_copy_number': scenario.g_copy_number,
                            'baseline_flux': scenario.baseline_flux,
                            'lethal_genes': len(scenario.top_genes[scenario.top_genes['impact'] == 'LETHAL']),
                            'lethal_rxns': len(scenario.top_rxns[scenario.top_rxns['impact'] == 'LETHAL']),
                            'significant_genes': len(scenario.top_genes[scenario.top_genes['impact'] == 'SIGNIFICANT']),
                            'significant_rxns': len(scenario.top_rxns[scenario.top_rxns['impact'] == 'SIGNIFICANT'])
                        })
                    else:
                        logger.warning(f"  Analysis incomplete for {scenario.scenario_id} (top_genes or top_rxns is None)")
                else:
                    logger.warning(f"  Baseline flux is zero or None for {scenario.scenario_id}, skipping summary")
                
            except Exception as e:
                logger.error(f"Error processing {scenario.scenario_id}: {e}")
                import traceback
                logger.error(traceback.format_exc())
                continue
        
        # Save scenario summary
        if scenario_summaries:
            summary_df = pd.DataFrame(scenario_summaries)
            summary_df.to_csv(output_dir / "scenario_summary.csv", index=False)
            logger.info(f"Scenario summary saved to: {output_dir / 'scenario_summary.csv'}")
        
        # =====================================================================
        # STEP 5: Compare scenarios and identify robust targets
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 5: COMPARING SCENARIOS")
        logger.info("=" * 70)
        
        comparison_results = compare_scenarios(scenarios)
        
        # Save robust targets
        comparison_results['robust_genes'].to_csv(
            output_dir / "robust_gene_targets.csv",
            index=False
        )
        comparison_results['robust_rxns'].to_csv(
            output_dir / "robust_reaction_targets.csv",
            index=False
        )
        
        logger.info(f"Robust targets saved:")
        logger.info(f"  Genes: {output_dir / 'robust_gene_targets.csv'}")
        logger.info(f"  Reactions: {output_dir / 'robust_reaction_targets.csv'}")
        
        # Generate comparison report
        generate_comparison_report(
            scenarios,
            comparison_results,
            output_dir / "comparison_report.txt"
        )
        
        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("SENSITIVITY ANALYSIS COMPLETE")
        logger.info("=" * 70)
        
        logger.info(f"\nResults Summary:")
        logger.info(f"  Total scenarios: {len(scenarios)}")
        logger.info(f"  Completed scenarios: {comparison_results['total_scenarios']}")
        logger.info(f"  Universal lethal genes: {comparison_results['universal_lethal_genes']}")
        logger.info(f"  Universal lethal reactions: {comparison_results['universal_lethal_rxns']}")
        logger.info(f"\nOutput directory: {output_dir}")
        
    except Exception as e:
        logger.error(f"Sensitivity analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

