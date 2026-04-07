"""
HMPV VBOF and Metabolic Model Construction Package
===================================================

This package provides modular tools for constructing a Viral Biomass Objective Function (VBOF)
for Human Metapneumovirus (HMPV) and integrating it with host genome-scale metabolic models.

Modules:
--------
- genome_analyzer: Parse and analyze HMPV genome composition
- protein_analyzer: Analyze viral protein sequences and amino acid composition
- energy_calculator: Calculate ATP/GTP costs for viral replication
- lipid_calculator: Calculate envelope lipid requirements
- vbof_builder: Construct the viral biomass objective function
- model_integration: Integrate VBOF with host metabolic models
- dual_objective_knockout: Dual-objective (BOF vs VBOF) gene knockout analysis
- analysis: Run FBA and identify antiviral targets

Author: Syed Mushahid Hussain
"""

__version__ = "1.1.0"
__author__ = "Syed Mushahid Hussain"

