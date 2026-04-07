"""
Dual-Objective Knockout Analysis Module
========================================

This module performs gene knockout analysis under two different objective functions:
1. Host biomass objective function (BOF) - maximizes host cell growth
2. Viral biomass objective function (VBOF) - maximizes virus production

The analysis compares the impact of each gene knockout on both objectives to
identify selective antiviral drug targets that hurt the virus while sparing the host.

Key Classes:
------------
- DualObjectiveConfig: Configuration for dual-objective analysis
- DualObjectiveKnockout: Main analysis class

Key Functions:
--------------
- perform_knockout_analysis: Run knockouts for a given objective
- merge_results: Combine host and virus knockout results
- classify_target: Classify genes based on impact thresholds
- filter_selective_targets: Find genes that selectively impact virus

Output Tables:
--------------
- host_growth_knockout_results.csv: Gene knockouts with BOF objective
- virus_growth_knockout_results.csv: Gene knockouts with VBOF objective
- merged_knockout_comparison.csv: Combined comparison table
- selective_antiviral_targets.csv: Filtered selective targets

Author: Syed Mushahid Hussain
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cobra
from cobra import Model
import pandas as pd
import numpy as np

from .config import (
    VBOF_REACTION_ID,
    HOST_BOF_REACTION_ID,
    DEFAULT_THRESHOLDS,
    ESSENTIALITY_THRESHOLD,
    SIGNIFICANT_THRESHOLD,
    MODERATE_THRESHOLD
)
from .exceptions import HMPVModelError

logger = logging.getLogger(__name__)


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class DualObjectiveConfig:
    """
    Configuration for dual-objective knockout analysis.
    
    Attributes:
    -----------
    vbof_id : str
        Reaction ID for viral biomass objective function (default: HMPV_VBOF)
    bof_id : str
        Reaction ID for host biomass objective function (default: R_biomass_hbec)
    thresholds : Dict
        Threshold configurations for target classification
    essentiality_threshold : float
        Flux ratio below which a knockout is considered lethal
    significant_threshold : float
        Flux ratio below which a knockout has significant impact
    moderate_threshold : float
        Flux ratio below which a knockout has moderate impact
    """
    vbof_id: str = VBOF_REACTION_ID
    bof_id: str = HOST_BOF_REACTION_ID
    thresholds: Dict[str, Dict[str, float]] = field(default_factory=lambda: DEFAULT_THRESHOLDS.copy())
    essentiality_threshold: float = ESSENTIALITY_THRESHOLD
    significant_threshold: float = SIGNIFICANT_THRESHOLD
    moderate_threshold: float = MODERATE_THRESHOLD


# =============================================================================
# TARGET CLASSIFICATION
# =============================================================================

def classify_impact(flux_ratio: float, config: DualObjectiveConfig) -> str:
    """
    Classify the impact of a knockout based on flux ratio.
    
    Parameters:
    -----------
    flux_ratio : float
        Ratio of knockout flux to baseline flux (0-1)
    config : DualObjectiveConfig
        Configuration with threshold values
    
    Returns:
    --------
    str : Impact classification (LETHAL, SIGNIFICANT, MODERATE, MINIMAL)
    """
    if flux_ratio < config.essentiality_threshold:
        return 'LETHAL'
    elif flux_ratio < config.significant_threshold:
        return 'SIGNIFICANT'
    elif flux_ratio < config.moderate_threshold:
        return 'MODERATE'
    else:
        return 'MINIMAL'


def classify_target(
    host_pct: float, 
    virus_pct: float, 
    thresholds: Dict[str, Dict[str, float]]
) -> str:
    """
    Classify gene knockout result into target categories.
    
    Parameters:
    -----------
    host_pct : float
        Host growth as percentage of wild-type (0-100)
    virus_pct : float
        Virus growth as percentage of wild-type (0-100)
    thresholds : Dict
        Threshold configurations for classification
    
    Returns:
    --------
    str : Target classification
        - CRITICAL_VIRAL_TARGET: Almost completely kills virus (may affect host)
        - HIGH_CONFIDENCE_TARGET: Strong effect on virus, minimal on host
        - SELECTIVE_TARGET: Good effect on virus, moderate on host
        - NON_SELECTIVE_TOXIC: Kills both virus and host
        - NOT_A_TARGET: Minimal effect on virus
    """
    # Convert percentages to fractions for comparison
    host_frac = host_pct / 100.0
    virus_frac = virus_pct / 100.0
    
    # Get threshold values
    lethal_virus = thresholds.get('lethal_virus', {})
    strict_selective = thresholds.get('strict_selective', {})
    selective = thresholds.get('selective_target', {})
    
    # Check for critical viral target (nearly zero virus growth)
    lethal_virus_max = lethal_virus.get('virus_max', 0.05)
    if virus_frac <= lethal_virus_max:
        return 'CRITICAL_VIRAL_TARGET'
    
    # Check for high confidence target (strict selective)
    strict_virus_max = strict_selective.get('virus_max', 0.1)
    strict_host_min = strict_selective.get('host_min', 0.9)
    if virus_frac <= strict_virus_max and host_frac >= strict_host_min:
        return 'HIGH_CONFIDENCE_TARGET'
    
    # Check for selective target
    selective_virus_max = selective.get('virus_max', 0.5)
    selective_host_min = selective.get('host_min', 0.8)
    if virus_frac <= selective_virus_max and host_frac >= selective_host_min:
        return 'SELECTIVE_TARGET'
    
    # Check for non-selective toxic (affects both)
    if virus_frac <= 0.5 and host_frac <= 0.5:
        return 'NON_SELECTIVE_TOXIC'
    
    return 'NOT_A_TARGET'


# =============================================================================
# DUAL-OBJECTIVE KNOCKOUT CLASS
# =============================================================================

class DualObjectiveKnockout:
    """
    Performs gene knockouts under both BOF and VBOF objectives.
    
    This class implements a comprehensive analysis pipeline that:
    1. Runs gene knockouts with host BOF as objective
    2. Runs gene knockouts with VBOF as objective
    3. Merges results into a comparison table
    4. Classifies targets based on selectivity thresholds
    5. Filters for selective antiviral targets
    
    Attributes:
    -----------
    model : cobra.Model
        The integrated metabolic model (host + VBOF)
    config : DualObjectiveConfig
        Configuration parameters
    host_baseline : float
        Baseline flux for host BOF (wild-type)
    virus_baseline : float
        Baseline flux for VBOF (wild-type)
    """
    
    def __init__(
        self, 
        model: Model, 
        config: Optional[DualObjectiveConfig] = None
    ):
        """
        Initialize the dual-objective knockout analyzer.
        
        Parameters:
        -----------
        model : cobra.Model
            The integrated metabolic model with both BOF and VBOF
        config : DualObjectiveConfig, optional
            Configuration parameters. Uses defaults if not provided.
        
        Raises:
        -------
        HMPVModelError : If required reactions are not found
        """
        self.model = model
        self.config = config or DualObjectiveConfig()
        
        # Validate model has required reactions
        self._validate_model()
        
        # Calculate baseline fluxes
        self.host_baseline = self._get_baseline_flux(self.config.bof_id)
        self.virus_baseline = self._get_baseline_flux(self.config.vbof_id)
        
        logger.info(f"Initialized DualObjectiveKnockout:")
        logger.info(f"  Host BOF ({self.config.bof_id}): baseline = {self.host_baseline:.6f}")
        logger.info(f"  VBOF ({self.config.vbof_id}): baseline = {self.virus_baseline:.6f}")
    
    def _validate_model(self) -> None:
        """
        Validate that the model has both required objective functions.
        
        Raises:
        -------
        HMPVModelError : If BOF or VBOF reaction not found
        """
        reaction_ids = [r.id for r in self.model.reactions]
        
        if self.config.bof_id not in reaction_ids:
            raise HMPVModelError(
                f"Host BOF reaction '{self.config.bof_id}' not found in model. "
                f"Available reactions containing 'biomass': "
                f"{[r for r in reaction_ids if 'biomass' in r.lower()]}"
            )
        
        if self.config.vbof_id not in reaction_ids:
            raise HMPVModelError(
                f"VBOF reaction '{self.config.vbof_id}' not found in model. "
                f"Please integrate HMPV VBOF first using integrate_model.py"
            )
        
        logger.info(f"Model validation passed: both BOF and VBOF found")
    
    def _get_baseline_flux(self, objective_id: str) -> float:
        """
        Calculate baseline flux for a given objective (wild-type).
        
        Parameters:
        -----------
        objective_id : str
            Reaction ID to use as objective
        
        Returns:
        --------
        float : Maximum flux (wild-type value)
        """
        with self.model:
            self.model.objective = objective_id
            solution = self.model.optimize()
            if solution.status == 'optimal':
                return solution.objective_value
            return 0.0
    
    def perform_knockout_analysis(
        self, 
        objective_id: str,
        baseline_flux: float,
        objective_name: str = "objective"
    ) -> pd.DataFrame:
        """
        Perform single gene knockout analysis for all genes with given objective.
        
        Parameters:
        -----------
        objective_id : str
            Reaction ID to use as objective function
        baseline_flux : float
            Wild-type flux for this objective
        objective_name : str
            Name for logging (e.g., "host" or "virus")
        
        Returns:
        --------
        pd.DataFrame : Knockout results with columns:
            - gene_id, gene_name, knockout_flux, baseline_flux
            - flux_ratio, growth_pct, impact, status
            - num_reactions, reactions
        """
        logger.info(f"Starting {objective_name} knockout analysis...")
        logger.info(f"  Objective: {objective_id}")
        logger.info(f"  Baseline flux: {baseline_flux:.6f}")
        logger.info(f"  Total genes: {len(self.model.genes)}")
        
        knockout_data = []
        total_genes = len(self.model.genes)
        
        for i, gene in enumerate(self.model.genes):
            if (i + 1) % 500 == 0:
                logger.info(f"  Progress: {i+1}/{total_genes} genes tested...")
            
            # Perform knockout using context manager
            with self.model:
                self.model.objective = objective_id
                gene.knock_out()
                try:
                    solution = self.model.optimize()
                    flux = solution.objective_value if solution.status == 'optimal' else 0.0
                    status = solution.status
                except Exception:
                    flux = 0.0
                    status = 'error'
            
            # Calculate metrics
            flux_ratio = flux / baseline_flux if baseline_flux > 0 else 0
            growth_pct = flux_ratio * 100
            
            # Classify impact
            impact = classify_impact(flux_ratio, self.config)
            
            # Get gene info
            associated_rxns = [r.id for r in gene.reactions]
            
            knockout_data.append({
                'gene_id': gene.id,
                'gene_name': gene.name if gene.name else gene.id,
                'knockout_flux': flux,
                'baseline_flux': baseline_flux,
                'flux_ratio': flux_ratio,
                'growth_pct': growth_pct,
                'impact': impact,
                'status': status,
                'num_reactions': len(associated_rxns),
                'reactions': ';'.join(associated_rxns[:10])
            })

        if not knockout_data:
            logger.warning(f"No knockout data collected for {objective_name} analysis.")
            return pd.DataFrame(columns=[
                'gene_id', 'gene_name', 'knockout_flux', 'baseline_flux',
                'flux_ratio', 'growth_pct', 'impact', 'status',
                'num_reactions', 'reactions'
            ])
        df = pd.DataFrame(knockout_data)
        df = df.sort_values('growth_pct', ascending=True)

        # Log summary
        logger.info(f"{objective_name.capitalize()} knockout analysis complete!")
        logger.info(f"  Lethal: {len(df[df['impact'] == 'LETHAL'])}")
        logger.info(f"  Significant: {len(df[df['impact'] == 'SIGNIFICANT'])}")
        logger.info(f"  Moderate: {len(df[df['impact'] == 'MODERATE'])}")
        logger.info(f"  Minimal: {len(df[df['impact'] == 'MINIMAL'])}")
        
        return df
    
    def run_host_knockout_analysis(self) -> pd.DataFrame:
        """
        Run gene knockout analysis with host BOF as objective.
        
        Returns:
        --------
        pd.DataFrame : Host knockout results
        """
        return self.perform_knockout_analysis(
            objective_id=self.config.bof_id,
            baseline_flux=self.host_baseline,
            objective_name="host"
        )
    
    def run_virus_knockout_analysis(self) -> pd.DataFrame:
        """
        Run gene knockout analysis with VBOF as objective.
        
        Returns:
        --------
        pd.DataFrame : Virus knockout results
        """
        return self.perform_knockout_analysis(
            objective_id=self.config.vbof_id,
            baseline_flux=self.virus_baseline,
            objective_name="virus"
        )
    
    def merge_results(
        self, 
        host_df: pd.DataFrame, 
        virus_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Merge host and virus knockout results into a comparison table.
        
        Parameters:
        -----------
        host_df : pd.DataFrame
            Host knockout results
        virus_df : pd.DataFrame
            Virus knockout results
        
        Returns:
        --------
        pd.DataFrame : Merged comparison table with columns:
            - gene_id, gene_name
            - host_flux, host_baseline, host_growth_pct, host_status
            - virus_flux, virus_baseline, virus_growth_pct, virus_status
            - selectivity_score, target_class
        """
        logger.info("Merging host and virus knockout results...")
        
        # Rename columns for merging
        host_renamed = host_df[['gene_id', 'gene_name', 'knockout_flux', 'baseline_flux', 
                                'growth_pct', 'impact', 'num_reactions', 'reactions']].copy()
        host_renamed.columns = ['gene_id', 'gene_name', 'host_flux', 'host_baseline',
                               'host_growth_pct', 'host_status', 'num_reactions', 'reactions']
        
        virus_renamed = virus_df[['gene_id', 'knockout_flux', 'baseline_flux', 
                                  'growth_pct', 'impact']].copy()
        virus_renamed.columns = ['gene_id', 'virus_flux', 'virus_baseline',
                                'virus_growth_pct', 'virus_status']
        
        # Merge on gene_id
        merged = pd.merge(host_renamed, virus_renamed, on='gene_id', how='outer')
        
        # Calculate selectivity score (higher = more selective for virus)
        # Positive score means host survives better than virus
        merged['selectivity_score'] = merged['host_growth_pct'] - merged['virus_growth_pct']
        
        # Classify each target
        merged['target_class'] = merged.apply(
            lambda row: classify_target(
                row['host_growth_pct'], 
                row['virus_growth_pct'],
                self.config.thresholds
            ),
            axis=1
        )
        
        # Sort by selectivity score (best targets first)
        merged = merged.sort_values('selectivity_score', ascending=False)
        
        logger.info(f"Merged results: {len(merged)} genes")
        logger.info(f"Target class distribution:")
        for tc, count in merged['target_class'].value_counts().items():
            logger.info(f"  {tc}: {count}")
        
        return merged
    
    def filter_selective_targets(
        self, 
        merged_df: pd.DataFrame,
        virus_max: Optional[float] = None,
        host_min: Optional[float] = None
    ) -> pd.DataFrame:
        """
        Filter for selective antiviral targets based on thresholds.
        
        Parameters:
        -----------
        merged_df : pd.DataFrame
            Merged knockout comparison table
        virus_max : float, optional
            Maximum virus growth percentage (default: 50%)
        host_min : float, optional
            Minimum host growth percentage (default: 80%)
        
        Returns:
        --------
        pd.DataFrame : Filtered selective targets
        """
        # Use default thresholds if not specified
        if virus_max is None:
            virus_max = self.config.thresholds.get('selective_target', {}).get('virus_max', 0.5) * 100
        if host_min is None:
            host_min = self.config.thresholds.get('selective_target', {}).get('host_min', 0.8) * 100
        
        logger.info(f"Filtering selective targets: virus < {virus_max}%, host > {host_min}%")
        
        # Apply filters
        selective = merged_df[
            (merged_df['virus_growth_pct'] <= virus_max) &
            (merged_df['host_growth_pct'] >= host_min)
        ].copy()
        
        # Sort by selectivity score
        selective = selective.sort_values('selectivity_score', ascending=False)
        
        logger.info(f"Found {len(selective)} selective targets")
        
        return selective
    
    def filter_critical_viral_targets(
        self, 
        merged_df: pd.DataFrame,
        virus_max: Optional[float] = None
    ) -> pd.DataFrame:
        """
        Filter for critical viral targets (nearly zero virus growth).
        
        These genes are essential for virus survival, even if host is affected.
        
        Parameters:
        -----------
        merged_df : pd.DataFrame
            Merged knockout comparison table
        virus_max : float, optional
            Maximum virus growth percentage (default: 5%)
        
        Returns:
        --------
        pd.DataFrame : Critical viral targets
        """
        if virus_max is None:
            virus_max = self.config.thresholds.get('lethal_virus', {}).get('virus_max', 0.05) * 100
        
        logger.info(f"Filtering critical viral targets: virus < {virus_max}%")
        
        critical = merged_df[merged_df['virus_growth_pct'] <= virus_max].copy()
        critical = critical.sort_values('host_growth_pct', ascending=False)
        
        logger.info(f"Found {len(critical)} critical viral targets")
        
        return critical
    
    def analyze_combined_objective(
        self,
        bof_weight: float = 0.9,
        vbof_weight: float = 0.1
    ) -> pd.DataFrame:
        """
        Analyze knockouts using weighted combined objective.
        
        This sets the objective to: bof_weight * BOF + vbof_weight * VBOF
        
        Parameters:
        -----------
        bof_weight : float
            Weight for host BOF (default: 0.9)
        vbof_weight : float
            Weight for VBOF (default: 0.1)
        
        Returns:
        --------
        pd.DataFrame : Combined objective knockout results
        """
        logger.info(f"Running combined objective analysis...")
        logger.info(f"  BOF weight: {bof_weight}")
        logger.info(f"  VBOF weight: {vbof_weight}")
        
        # Calculate combined baseline
        combined_baseline = (
            bof_weight * self.host_baseline + 
            vbof_weight * self.virus_baseline
        )
        
        knockout_data = []
        total_genes = len(self.model.genes)
        
        for i, gene in enumerate(self.model.genes):
            if (i + 1) % 500 == 0:
                logger.info(f"  Progress: {i+1}/{total_genes} genes tested...")
            
            with self.model:
                # Set combined objective
                self.model.objective = {
                    self.model.reactions.get_by_id(self.config.bof_id): bof_weight,
                    self.model.reactions.get_by_id(self.config.vbof_id): vbof_weight
                }
                
                gene.knock_out()
                try:
                    solution = self.model.optimize()
                    if solution.status == 'optimal':
                        # Get individual fluxes
                        bof_flux = solution.fluxes.get(self.config.bof_id, 0)
                        vbof_flux = solution.fluxes.get(self.config.vbof_id, 0)
                        combined_flux = bof_weight * bof_flux + vbof_weight * vbof_flux
                        status = 'optimal'
                    else:
                        bof_flux = vbof_flux = combined_flux = 0.0
                        status = solution.status
                except Exception:
                    bof_flux = vbof_flux = combined_flux = 0.0
                    status = 'error'
            
            # Calculate metrics
            combined_ratio = combined_flux / combined_baseline if combined_baseline > 0 else 0
            host_ratio = bof_flux / self.host_baseline if self.host_baseline > 0 else 0
            virus_ratio = vbof_flux / self.virus_baseline if self.virus_baseline > 0 else 0
            
            knockout_data.append({
                'gene_id': gene.id,
                'gene_name': gene.name if gene.name else gene.id,
                'host_flux': bof_flux,
                'virus_flux': vbof_flux,
                'combined_flux': combined_flux,
                'combined_baseline': combined_baseline,
                'host_growth_pct': host_ratio * 100,
                'virus_growth_pct': virus_ratio * 100,
                'combined_growth_pct': combined_ratio * 100,
                'bof_weight': bof_weight,
                'vbof_weight': vbof_weight,
                'status': status
            })

        if not knockout_data:
            logger.warning("No knockout data collected for combined objective analysis.")
            return pd.DataFrame(columns=[
                'gene_id', 'gene_name', 'bof_flux', 'vbof_flux', 'bof_ratio',
                'vbof_ratio', 'combined_ratio', 'combined_growth_pct',
                'bof_weight', 'vbof_weight', 'status'
            ])
        df = pd.DataFrame(knockout_data)
        df = df.sort_values('combined_growth_pct', ascending=True)

        logger.info(f"Combined objective analysis complete!")
        
        return df
    
    def run_full_analysis(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run the complete dual-objective knockout analysis pipeline.
        
        Returns:
        --------
        Tuple of DataFrames:
            - host_results: Host knockout results
            - virus_results: Virus knockout results
            - merged_results: Merged comparison table
            - selective_targets: Selective antiviral targets
            - critical_targets: Critical viral targets
        """
        logger.info("=" * 70)
        logger.info("RUNNING FULL DUAL-OBJECTIVE ANALYSIS")
        logger.info("=" * 70)
        
        # Step 1: Run host knockout analysis
        logger.info("\n--- Step 1: Host Knockout Analysis ---")
        host_results = self.run_host_knockout_analysis()
        
        # Step 2: Run virus knockout analysis
        logger.info("\n--- Step 2: Virus Knockout Analysis ---")
        virus_results = self.run_virus_knockout_analysis()
        
        # Step 3: Merge results
        logger.info("\n--- Step 3: Merging Results ---")
        merged_results = self.merge_results(host_results, virus_results)
        
        # Step 4: Filter selective targets
        logger.info("\n--- Step 4: Filtering Selective Targets ---")
        selective_targets = self.filter_selective_targets(merged_results)
        
        # Step 5: Filter critical viral targets
        logger.info("\n--- Step 5: Filtering Critical Viral Targets ---")
        critical_targets = self.filter_critical_viral_targets(merged_results)
        
        logger.info("\n" + "=" * 70)
        logger.info("ANALYSIS COMPLETE")
        logger.info("=" * 70)
        
        return host_results, virus_results, merged_results, selective_targets, critical_targets

    # =========================================================================
    # REACTION KNOCKOUT ANALYSIS
    # =========================================================================
    
    def perform_reaction_knockout_analysis(
        self,
        objective_id: str,
        baseline_flux: float,
        objective_name: str = "objective"
    ) -> pd.DataFrame:
        """
        Perform single reaction knockout analysis for all reactions with given objective.
        
        Parameters:
        -----------
        objective_id : str
            Reaction ID to use as objective function
        baseline_flux : float
            Wild-type flux for this objective
        objective_name : str
            Name for logging (e.g., "host" or "virus")
        
        Returns:
        --------
        pd.DataFrame : Reaction knockout results with columns:
            - reaction_id, reaction_name, subsystem
            - knockout_flux, baseline_flux, flux_ratio, growth_pct
            - impact, status, num_genes, gene_ids, gene_names, equation
        """
        logger.info(f"Starting {objective_name} REACTION knockout analysis...")
        logger.info(f"  Objective: {objective_id}")
        logger.info(f"  Baseline flux: {baseline_flux:.6f}")
        
        # Exclude exchange reactions and objective reactions
        reactions_to_test = [
            r for r in self.model.reactions 
            if not r.id.startswith('EX_') 
            and r.id != self.config.vbof_id
            and r.id != self.config.bof_id
        ]
        
        logger.info(f"  Total reactions to test: {len(reactions_to_test)}")
        
        knockout_data = []
        total_rxns = len(reactions_to_test)
        
        for i, rxn in enumerate(reactions_to_test):
            if (i + 1) % 500 == 0:
                logger.info(f"  Progress: {i+1}/{total_rxns} reactions tested...")
            
            # Perform knockout using context manager
            with self.model:
                self.model.objective = objective_id
                rxn.knock_out()
                try:
                    solution = self.model.optimize()
                    flux = solution.objective_value if solution.status == 'optimal' else 0.0
                    status = solution.status
                except Exception:
                    flux = 0.0
                    status = 'error'
            
            # Calculate metrics
            flux_ratio = flux / baseline_flux if baseline_flux > 0 else 0
            growth_pct = flux_ratio * 100
            
            # Classify impact
            impact = classify_impact(flux_ratio, self.config)
            
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
                'growth_pct': growth_pct,
                'impact': impact,
                'status': status,
                'num_genes': len(associated_genes),
                'gene_ids': ';'.join(associated_genes[:10]),
                'gene_names': ';'.join(gene_names[:10]),
                'equation': rxn.reaction
            })

        if not knockout_data:
            logger.warning(f"No reaction knockout data collected for {objective_name} analysis.")
            return pd.DataFrame(columns=[
                'reaction_id', 'reaction_name', 'knockout_flux', 'baseline_flux',
                'flux_ratio', 'growth_pct', 'impact', 'status',
                'num_genes', 'gene_ids', 'gene_names', 'equation'
            ])
        df = pd.DataFrame(knockout_data)
        df = df.sort_values('growth_pct', ascending=True)

        # Log summary
        logger.info(f"{objective_name.capitalize()} REACTION knockout analysis complete!")
        logger.info(f"  Lethal: {len(df[df['impact'] == 'LETHAL'])}")
        logger.info(f"  Significant: {len(df[df['impact'] == 'SIGNIFICANT'])}")
        logger.info(f"  Moderate: {len(df[df['impact'] == 'MODERATE'])}")
        logger.info(f"  Minimal: {len(df[df['impact'] == 'MINIMAL'])}")
        
        return df
    
    def run_host_reaction_knockout_analysis(self) -> pd.DataFrame:
        """
        Run reaction knockout analysis with host BOF as objective.
        
        Returns:
        --------
        pd.DataFrame : Host reaction knockout results
        """
        return self.perform_reaction_knockout_analysis(
            objective_id=self.config.bof_id,
            baseline_flux=self.host_baseline,
            objective_name="host"
        )
    
    def run_virus_reaction_knockout_analysis(self) -> pd.DataFrame:
        """
        Run reaction knockout analysis with VBOF as objective.
        
        Returns:
        --------
        pd.DataFrame : Virus reaction knockout results
        """
        return self.perform_reaction_knockout_analysis(
            objective_id=self.config.vbof_id,
            baseline_flux=self.virus_baseline,
            objective_name="virus"
        )
    
    def merge_reaction_results(
        self,
        host_df: pd.DataFrame,
        virus_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Merge host and virus reaction knockout results into a comparison table.
        
        Parameters:
        -----------
        host_df : pd.DataFrame
            Host reaction knockout results
        virus_df : pd.DataFrame
            Virus reaction knockout results
        
        Returns:
        --------
        pd.DataFrame : Merged comparison table with columns:
            - reaction_id, reaction_name, subsystem
            - host_flux, host_baseline, host_growth_pct, host_status
            - virus_flux, virus_baseline, virus_growth_pct, virus_status
            - selectivity_score, target_class
        """
        logger.info("Merging host and virus REACTION knockout results...")
        
        # Rename columns for merging
        host_renamed = host_df[['reaction_id', 'reaction_name', 'subsystem', 'knockout_flux', 
                                'baseline_flux', 'growth_pct', 'impact', 'num_genes', 
                                'gene_ids', 'gene_names', 'equation']].copy()
        host_renamed.columns = ['reaction_id', 'reaction_name', 'subsystem', 'host_flux', 
                               'host_baseline', 'host_growth_pct', 'host_status', 'num_genes',
                               'gene_ids', 'gene_names', 'equation']
        
        virus_renamed = virus_df[['reaction_id', 'knockout_flux', 'baseline_flux', 
                                  'growth_pct', 'impact']].copy()
        virus_renamed.columns = ['reaction_id', 'virus_flux', 'virus_baseline',
                                'virus_growth_pct', 'virus_status']
        
        # Merge on reaction_id
        merged = pd.merge(host_renamed, virus_renamed, on='reaction_id', how='outer')
        
        # Calculate selectivity score (higher = more selective for virus)
        merged['selectivity_score'] = merged['host_growth_pct'] - merged['virus_growth_pct']
        
        # Classify each target
        merged['target_class'] = merged.apply(
            lambda row: classify_target(
                row['host_growth_pct'], 
                row['virus_growth_pct'],
                self.config.thresholds
            ),
            axis=1
        )
        
        # Sort by selectivity score (best targets first)
        merged = merged.sort_values('selectivity_score', ascending=False)
        
        logger.info(f"Merged reaction results: {len(merged)} reactions")
        logger.info(f"Target class distribution:")
        for tc, count in merged['target_class'].value_counts().items():
            logger.info(f"  {tc}: {count}")
        
        return merged
    
    def filter_selective_reaction_targets(
        self,
        merged_df: pd.DataFrame,
        virus_max: Optional[float] = None,
        host_min: Optional[float] = None
    ) -> pd.DataFrame:
        """
        Filter for selective antiviral reaction targets based on thresholds.
        
        Parameters:
        -----------
        merged_df : pd.DataFrame
            Merged reaction knockout comparison table
        virus_max : float, optional
            Maximum virus growth percentage (default: 50%)
        host_min : float, optional
            Minimum host growth percentage (default: 80%)
        
        Returns:
        --------
        pd.DataFrame : Filtered selective reaction targets
        """
        if virus_max is None:
            virus_max = self.config.thresholds.get('selective_target', {}).get('virus_max', 0.5) * 100
        if host_min is None:
            host_min = self.config.thresholds.get('selective_target', {}).get('host_min', 0.8) * 100
        
        logger.info(f"Filtering selective REACTION targets: virus < {virus_max}%, host > {host_min}%")
        
        selective = merged_df[
            (merged_df['virus_growth_pct'] <= virus_max) &
            (merged_df['host_growth_pct'] >= host_min)
        ].copy()
        
        selective = selective.sort_values('selectivity_score', ascending=False)
        
        logger.info(f"Found {len(selective)} selective reaction targets")
        
        return selective
    
    def filter_critical_reaction_targets(
        self,
        merged_df: pd.DataFrame,
        virus_max: Optional[float] = None
    ) -> pd.DataFrame:
        """
        Filter for critical viral reaction targets (nearly zero virus growth).
        
        Parameters:
        -----------
        merged_df : pd.DataFrame
            Merged reaction knockout comparison table
        virus_max : float, optional
            Maximum virus growth percentage (default: 5%)
        
        Returns:
        --------
        pd.DataFrame : Critical viral reaction targets
        """
        if virus_max is None:
            virus_max = self.config.thresholds.get('lethal_virus', {}).get('virus_max', 0.05) * 100
        
        logger.info(f"Filtering critical viral REACTION targets: virus < {virus_max}%")
        
        critical = merged_df[merged_df['virus_growth_pct'] <= virus_max].copy()
        critical = critical.sort_values('host_growth_pct', ascending=False)
        
        logger.info(f"Found {len(critical)} critical viral reaction targets")
        
        return critical
    
    def analyze_subsystem_essentiality(
        self,
        merged_reaction_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Analyze which metabolic subsystems contain selective antiviral targets.
        
        Parameters:
        -----------
        merged_reaction_df : pd.DataFrame
            Merged reaction knockout comparison table
        
        Returns:
        --------
        pd.DataFrame : Subsystem essentiality analysis with dual-objective metrics
        """
        logger.info("Analyzing subsystem essentiality for dual-objective analysis...")
        
        subsystem_data = []
        
        for subsystem in merged_reaction_df['subsystem'].unique():
            subsystem_rxns = merged_reaction_df[merged_reaction_df['subsystem'] == subsystem]
            
            # Count by target class
            critical_count = len(subsystem_rxns[subsystem_rxns['target_class'] == 'CRITICAL_VIRAL_TARGET'])
            high_conf_count = len(subsystem_rxns[subsystem_rxns['target_class'] == 'HIGH_CONFIDENCE_TARGET'])
            selective_count = len(subsystem_rxns[subsystem_rxns['target_class'] == 'SELECTIVE_TARGET'])
            non_selective_count = len(subsystem_rxns[subsystem_rxns['target_class'] == 'NON_SELECTIVE_TOXIC'])
            not_target_count = len(subsystem_rxns[subsystem_rxns['target_class'] == 'NOT_A_TARGET'])
            total_count = len(subsystem_rxns)
            
            # Calculate averages
            avg_host_growth = subsystem_rxns['host_growth_pct'].mean()
            avg_virus_growth = subsystem_rxns['virus_growth_pct'].mean()
            avg_selectivity = subsystem_rxns['selectivity_score'].mean()
            max_selectivity = subsystem_rxns['selectivity_score'].max()
            
            # Get sample selective reactions
            selective_rxns = subsystem_rxns[
                subsystem_rxns['target_class'].isin(['SELECTIVE_TARGET', 'HIGH_CONFIDENCE_TARGET', 'CRITICAL_VIRAL_TARGET'])
            ]['reaction_name'].head(3).tolist()
            
            # Calculate essentiality score (prioritizes selective targets)
            essentiality_score = (
                high_conf_count * 3 + 
                selective_count * 2 + 
                critical_count * 1
            ) / total_count if total_count > 0 else 0
            
            subsystem_data.append({
                'subsystem': subsystem,
                'total_reactions': total_count,
                'critical_viral_targets': critical_count,
                'high_confidence_targets': high_conf_count,
                'selective_targets': selective_count,
                'non_selective_toxic': non_selective_count,
                'not_a_target': not_target_count,
                'avg_host_growth_pct': avg_host_growth,
                'avg_virus_growth_pct': avg_virus_growth,
                'avg_selectivity_score': avg_selectivity,
                'max_selectivity_score': max_selectivity,
                'essentiality_score': essentiality_score,
                'sample_selective_reactions': ';'.join(selective_rxns)
            })
        
        df = pd.DataFrame(subsystem_data)
        df = df.sort_values('essentiality_score', ascending=False)
        
        logger.info(f"Analyzed {len(df)} subsystems")
        
        return df
    
    def run_full_reaction_analysis(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run the complete dual-objective REACTION knockout analysis pipeline.
        
        Returns:
        --------
        Tuple of DataFrames:
            - host_results: Host reaction knockout results
            - virus_results: Virus reaction knockout results
            - merged_results: Merged comparison table
            - selective_targets: Selective reaction targets
            - critical_targets: Critical viral reaction targets
            - subsystem_analysis: Subsystem essentiality analysis
        """
        logger.info("=" * 70)
        logger.info("RUNNING FULL DUAL-OBJECTIVE REACTION ANALYSIS")
        logger.info("=" * 70)
        
        # Step 1: Run host reaction knockout analysis
        logger.info("\n--- Step 1: Host Reaction Knockout Analysis ---")
        host_results = self.run_host_reaction_knockout_analysis()
        
        # Step 2: Run virus reaction knockout analysis
        logger.info("\n--- Step 2: Virus Reaction Knockout Analysis ---")
        virus_results = self.run_virus_reaction_knockout_analysis()
        
        # Step 3: Merge results
        logger.info("\n--- Step 3: Merging Reaction Results ---")
        merged_results = self.merge_reaction_results(host_results, virus_results)
        
        # Step 4: Filter selective targets
        logger.info("\n--- Step 4: Filtering Selective Reaction Targets ---")
        selective_targets = self.filter_selective_reaction_targets(merged_results)
        
        # Step 5: Filter critical viral targets
        logger.info("\n--- Step 5: Filtering Critical Viral Reaction Targets ---")
        critical_targets = self.filter_critical_reaction_targets(merged_results)
        
        # Step 6: Subsystem analysis
        logger.info("\n--- Step 6: Subsystem Essentiality Analysis ---")
        subsystem_analysis = self.analyze_subsystem_essentiality(merged_results)
        
        logger.info("\n" + "=" * 70)
        logger.info("REACTION ANALYSIS COMPLETE")
        logger.info("=" * 70)
        
        return host_results, virus_results, merged_results, selective_targets, critical_targets, subsystem_analysis


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_integrated_model(model_path: Path, vbof_id: str = VBOF_REACTION_ID) -> Model:
    """
    Load the integrated model and validate it has VBOF.
    
    Parameters:
    -----------
    model_path : Path
        Path to the integrated model file
    vbof_id : str
        Expected VBOF reaction ID
    
    Returns:
    --------
    cobra.Model : Loaded model
    
    Raises:
    -------
    FileNotFoundError : If model file not found
    HMPVModelError : If VBOF not found in model
    """
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    logger.info(f"Loading model from: {model_path}")
    
    # Suppress COBRApy warnings
    cobra_logger = logging.getLogger('cobra.io.sbml')
    original_level = cobra_logger.level
    cobra_logger.setLevel(logging.CRITICAL)
    
    try:
        model = cobra.io.read_sbml_model(str(model_path))
    finally:
        cobra_logger.setLevel(original_level)
    
    logger.info(f"Model loaded: {len(model.reactions)} reactions, {len(model.genes)} genes")

    return model
