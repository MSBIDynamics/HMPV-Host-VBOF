[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Repo Size](https://img.shields.io/github/repo-size/MSBIDynamics/HMPV-Host-VBOF)
![GitHub release (with filter)](https://img.shields.io/github/v/release/MSBIDynamics/HMPV-Host-VBOF?logo=github&label=HMPV-Host-VBOFs&color=B4A069&style=flat-square)
[![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8270303-B4A069?style=flat-square&logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.20628356)

# HMPV Metabolic Modeling Project

## Table of Contents

<img align="right" src="https://github.com/MSBIDynamics/HMPV-Host-VBOF/raw/main/Pipeline_HMPV_VBOF.png" height="400" width="400"/>

-----
1. [Project Overview](#project-overview)
2. [Key Research Findings](#key-research-findings)
3. [Project Structure](#project-structure)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Workflow Pipeline](#workflow-pipeline)
7. [Module Documentation](#module-documentation)
8. [Output Files](#output-files)
9. [Sensitivity Analysis](#sensitivity-analysis)
10. [Configuration](#configuration)
11. [References](#references)

---

## Project Overview

This project constructs a **Viral Biomass Objective Function (VBOF)** for **Human Metapneumovirus (HMPV)** and integrates it with host genome-scale metabolic models (GEMs) to identify potential antiviral drug targets.

### Key Features

- **Modular VBOF Construction**: Builds VBOF from genome and protein sequences
- **Host Model Integration**: Integrates VBOF with human metabolic models (HBEC)
- **Dual-Objective Knockout Analysis**: Novel approach evaluating both host viability AND viral production
- **Selective Target Identification**: Identifies targets that kill virus while sparing host cells
- **Sensitivity Analysis**: Tests robustness of targets across different VBOF parameter combinations
- **Comprehensive Documentation**: Well-documented code following PEP8 standards

### Scientific Background

The VBOF represents the stoichiometric requirements for producing one HMPV virion, including:
- **Genome**: Nucleotide requirements (ATP, UTP, GTP, CTP)
- **Proteins**: Amino acid requirements based on protein copy numbers
- **Lipids**: Envelope lipid composition (phospholipids, cholesterol, sphingomyelin)
- **Glycans**: N-linked and O-linked glycosylation requirements
- **Energy**: ATP/GTP costs for replication and translation

---

## Key Research Findings

### Primary Discovery: Two Highly Selective Antiviral Targets

Our dual-objective knockout analysis of 1,390 host genes identified **two ideal drug targets**:

| Gene | Protein | Host Growth | Virus Growth | Selectivity Score |
|------|---------|-------------|--------------|-------------------|
| **PGM3** | Phosphoglucomutase 3 | 100% | 0% | **100.0** (perfect) |
| **GNPNAT1** | Glucosamine-phosphate N-acetyltransferase 1 | 100% | 2.6% | **97.4** |

### Target Pathway: Hexosamine Biosynthesis

Both targets catalyze sequential reactions in the **hexosamine biosynthesis pathway (HBP)**, which produces UDP-N-acetylglucosamine (UDP-GlcNAc) - the essential donor for glycosylation of viral F and G surface proteins.

```
Glucose-6-P → Fructose-6-P → Glucosamine-6-P
                                    ↓ GNPNAT1 ★
                              N-Acetylglucosamine-6-P
                                    ↓ PGM3 ★
                              N-Acetylglucosamine-1-P
                                    ↓
                              UDP-N-Acetylglucosamine
                                    ↓
                              Viral Glycoprotein Glycosylation
```

### Why This Works

- **HMPV requires ~13,300 glycosylation events per virion** (F and G proteins)
- **Host cells tolerate HBP disruption** due to lower steady-state demands and glycan recycling
- **Viral production creates acute metabolic bottleneck** at HBP that exceeds host capacity

### Complete Analysis Summary

| Metric | Value |
|--------|-------|
| Genes tested | 1,390 |
| Reactions tested | 1,769 |
| Selective gene targets | 2 |
| Selective reaction targets | 7 |
| Critical viral targets | 29 genes, 58 reactions |
| VBOF baseline flux | 8.68 virions/cell/h |

For detailed findings, see [RESEARCH_PAPER.md](RESEARCH_PAPER.md) and [RESEARCH_HIGHLIGHTS.md](RESEARCH_HIGHLIGHTS.md)

---

## Project Structure

```
GEM_HMPV-Host_Interactions/
├── src/                          # Core modules
│   ├── __init__.py              # Package initialization
│   ├── config.py                 # Configuration parameters and constants
│   ├── exceptions.py             # Custom exception classes
│   ├── genome_analyzer.py        # Genome sequence analysis
│   ├── protein_analyzer.py       # Protein sequence analysis
│   ├── glycan_analyzer.py        # Glycosylation analysis
│   ├── vbof_builder.py           # VBOF construction
│   └── model_integration.py       # Model integration with COBRApy
│
├── Data/                         # Input data directory
│   ├── genomic/                  # Genome sequence files
│   │   ├── GCF_002815375.1_ASM281537v1_genomic.fna
│   │   └── GCF_002815375.1_ASM281537v1_genomic.gff
│   ├── protein/                  # Protein sequence files
│   │   └── GCF_002815375.1_ASM281537v1_protein.faa
│   └── sbml/                     # Host metabolic models (SBML format)
│       ├── iHsaEC21.xml
│       └── iHsaEC21_clean.xml
│
├── docs/                         # Documentation
│   ├── research_knowledge_base.md  # Literature review and methodology
│   └── quick_reference.md          # Quick reference guide
│
├── output/                       # Generated output files
│   ├── hmpv_vbof.json            # Raw VBOF stoichiometry
│   ├── hmpv_vbof_normalized.json # Normalized VBOF
│   ├── hmpv_vbof_summary.txt     # Human-readable VBOF summary
│   ├── iHsaEC21_CLEAN_with_HMPV_VBOF.xml  # Integrated model
│   ├── antiviral_analysis/       # Antiviral target analysis results
│   └── sensitivity_analysis/     # Sensitivity analysis results
│
├── build_vbof.py                 # Main VBOF construction script
├── normalize_vbof.py             # VBOF normalization script
├── integrate_model.py            # Model integration script
├── clean_host_model.py            # Host model cleaning script
├── antiviral_target_analysis.py  # Antiviral target identification
├── vbof_sensitivity_analysis.py  # Sensitivity analysis across parameters
├── generate_report.py            # Report generation
├── requirements.txt              # Python dependencies
└── README.md                     # This file
```

---

## Installation

### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

### Step 1: Clone or Download the Project

```bash
cd /path/to/GEM_HMPV-Host_Interactions
```

### Step 2: Install Dependencies

```bash
pip install -r requirements.txt
```

### Required Packages

- **numpy** (≥1.21.0): Numerical computations
- **pandas** (≥1.3.0): Data manipulation and analysis
- **scipy** (≥1.7.0): Scientific computing
- **cobra** (≥0.26.0): Constraint-based metabolic modeling (COBRApy)

### Step 3: Verify Installation

```bash
python -c "import cobra; import pandas; import numpy; print('All packages installed successfully!')"
```

---

## Quick Start

### Basic Workflow

1. **Build VBOF** from genome and protein data
2. **Normalize VBOF** coefficients for FBA
3. **Integrate VBOF** into host metabolic model
4. **Run Antiviral Analysis** to identify targets
5. **Perform Sensitivity Analysis** to test robustness

```bash
# Step 1: Build VBOF
python build_vbof.py

# Step 2: Normalize VBOF
python normalize_vbof.py

# Step 3: Integrate into host model
python integrate_model.py

# Step 4: Identify antiviral targets
python antiviral_target_analysis.py

# Step 5: Sensitivity analysis (optional)
python vbof_sensitivity_analysis.py
```

---

## Workflow Pipeline

### 1. VBOF Construction (`build_vbof.py`)

**Purpose**: Constructs the Viral Biomass Objective Function from HMPV genome and protein sequences.

**What it does**:
- Loads HMPV genome FASTA file
- Parses GFF annotations
- Calculates nucleotide stoichiometry
- Loads protein sequences
- Calculates amino acid stoichiometry with copy numbers
- Builds complete VBOF including:
  - Genome nucleotides
  - Protein amino acids
  - Lipid envelope components
  - Glycan requirements
  - Energy costs (ATP/GTP)

**Output**:
- `output/hmpv_vbof.json`: Complete VBOF stoichiometry (raw coefficients)
- `output/hmpv_vbof_summary.txt`: Human-readable summary

**Usage**:
```bash
python build_vbof.py
```

**Key Parameters** (in `src/config.py`):
- `VIRION_DIAMETER_NM`: Virion diameter in nanometers (default: 200.0)
- `HMPV_COPY_NUMBERS`: Protein copy numbers per virion
- `LIPID_COMPOSITION`: Envelope lipid fractions

---

### 2. VBOF Normalization (`normalize_vbof.py`)

**Purpose**: Normalizes VBOF coefficients to reasonable values for flux balance analysis.

**What it does**:
- Loads raw VBOF from `hmpv_vbof.json`
- Normalizes coefficients by dividing by total consumed
- Ensures coefficients are in range suitable for FBA (-1 to +1 for major components)

**Why needed**: Raw VBOF has coefficients in millions (molecule counts), which makes FBA infeasible.

**Output**:
- `output/hmpv_vbof_normalized.json`: Normalized VBOF stoichiometry

**Usage**:
```bash
python normalize_vbof.py
```

---

### 3. Model Integration (`integrate_model.py`)

**Purpose**: Integrates the normalized VBOF into the host metabolic model.

**What it does**:
- Loads normalized VBOF
- Loads host metabolic model (iHsaEC21)
- Maps VBOF metabolite IDs to host model IDs
- Creates VBOF reaction in COBRApy format
- Adds reaction to model
- Validates integrated model
- Performs test FBA

**Output**:
- `output/iHsaEC21_CLEAN_with_HMPV_VBOF.xml`: Integrated model (SBML format)
- `output/integration_summary.txt`: Integration summary and validation results

**Usage**:
```bash
python integrate_model.py
```

**Note**: The integrated model sets VBOF as the objective function for FBA.

---

### 4. Antiviral Target Analysis (`antiviral_target_analysis.py`)

**Purpose**: Identifies host genes and reactions essential for HMPV replication.

**What it does**:
- Loads integrated model
- Calculates baseline VBOF flux (wild-type)
- Performs single gene knockout analysis
- Performs single reaction knockout analysis
- Classifies impact levels:
  - **LETHAL**: >99% flux reduction
  - **SIGNIFICANT**: 50-99% reduction
  - **MODERATE**: 10-50% reduction
  - **MINIMAL**: <10% reduction
- Analyzes subsystem essentiality
- Generates comprehensive reports

**Output** (in `output/antiviral_analysis/`):
- `gene_knockout_results.csv`: Complete gene knockout results
- `reaction_knockout_results.csv`: Complete reaction knockout results
- `top_gene_targets.csv`: Top antiviral gene targets
- `top_reaction_targets.csv`: Top antiviral reaction targets
- `subsystem_essentiality.csv`: Subsystem-level analysis
- `antiviral_targets_report.txt`: Comprehensive text report

**Usage**:
```bash
python antiviral_target_analysis.py
```

**Runtime**: This analysis can take 30-60 minutes depending on model size.

---

### 5. Sensitivity Analysis (`vbof_sensitivity_analysis.py`)

**Purpose**: Tests robustness of antiviral targets across different VBOF parameter combinations.

**What it does**:
- Generates multiple VBOF variants with different parameters:
  - Virion diameter: 150, 200, 250, 300 nm
  - F protein copy number: 250, 350, 450
  - G protein copy number: 150, 250, 350
- Builds and integrates VBOF for each scenario (36 total)
- Runs antiviral analysis for each scenario
- Compares results across scenarios
- Identifies robust targets (appear in multiple scenarios)
- Calculates robustness scores

**Output** (in `output/sensitivity_analysis/`):
- `scenario_XXX/`: Individual results for each scenario
- `scenario_summary.csv`: Summary of all scenarios
- `robust_gene_targets.csv`: Genes ranked by robustness
- `robust_reaction_targets.csv`: Reactions ranked by robustness
- `comparison_report.txt`: Cross-scenario comparison report

**Usage**:
```bash
python vbof_sensitivity_analysis.py
```

**Runtime**: This analysis can take several hours (36 scenarios × ~30-60 min each).

**Interpretation**:
- **Universal targets** (100% robustness): Appear in ALL scenarios - highest priority
- **High robustness** (>70%): Reliable across most parameter combinations
- **Moderate robustness** (30-70%): Sensitive to parameter uncertainty
- **Low robustness** (<30%): May be scenario-specific

---

### 6. Dual-Objective Analysis (`dual_objective_analysis.py`)

**Purpose**: Identifies selective antiviral drug targets by comparing host growth (BOF) vs virus production (VBOF) for each gene knockout. This is the **key methodological innovation** of this project.

**Why Dual-Objective?**
Previous metabolic modeling studies evaluated only viral production knockouts. This approach cannot distinguish:
- Targets that kill both virus AND host (toxic)
- Targets that kill virus but spare host (selective - ideal drug targets)

Our dual-objective framework solves this by evaluating each knockout under BOTH objectives.

**What it does**:
- Loads integrated model (host + HMPV VBOF)
- Runs gene knockouts with host BOF as objective → measures impact on host growth
- Runs gene knockouts with VBOF as objective → measures impact on virus production
- Runs reaction knockouts with both objectives
- Calculates selectivity scores: `Selectivity = Host Growth % - Virus Growth %`
- Classifies targets into five categories
- Generates comprehensive reports with pathway analysis

**Output** (in `output_HBEC/dual_objective_analysis/`):

*Gene-level analysis:*
- `host_growth_knockout_results.csv`: All 1,390 gene knockouts with host BOF objective
- `virus_growth_knockout_results.csv`: All gene knockouts with VBOF objective
- `merged_knockout_comparison.csv`: Combined comparison table
- `selective_antiviral_targets.csv`: **Key result**: 2 selective targets (PGM3, GNPNAT1)
- `critical_viral_targets.csv`: 29 critical viral targets

*Reaction-level analysis:*
- `host_reaction_knockout_results.csv`: All 1,769 reaction knockouts (host objective)
- `virus_reaction_knockout_results.csv`: Reaction knockouts (virus objective)
- `merged_reaction_knockout_comparison.csv`: Combined reaction comparison
- `selective_reaction_targets.csv`: **Key result**: 7 selective reaction targets
- `critical_reaction_targets.csv`: 58 critical viral reactions
- `reaction_subsystem_essentiality.csv`: Pathway-level analysis

*Report:*
- `dual_objective_report.txt`: Comprehensive human-readable report

**Usage**:
```bash
# Default analysis (virus <50%, host >80%)
python dual_objective_analysis.py

# Strict thresholds (virus <10%, host >90%)
python dual_objective_analysis.py --virus-max 10 --host-min 90

# With combined objective analysis
python dual_objective_analysis.py --combined --bof-weight 0.9 --vbof-weight 0.1
```

**Command-Line Options**:
- `--virus-max FLOAT`: Maximum virus growth % for selective targets (default: 50)
- `--host-min FLOAT`: Minimum host growth % for selective targets (default: 80)
- `--critical FLOAT`: Maximum virus growth % for critical targets (default: 5)
- `--combined`: Run combined objective analysis
- `--bof-weight FLOAT`: BOF weight for combined analysis (default: 0.9)
- `--vbof-weight FLOAT`: VBOF weight for combined analysis (default: 0.1)

**Target Classification**:

| Classification | Virus Growth | Host Growth | Interpretation |
|----------------|--------------|-------------|----------------|
| **HIGH_CONFIDENCE_TARGET** | <10% | >90% | Ideal drug targets |
| **SELECTIVE_TARGET** | <50% | >80% | Good drug candidates |
| **CRITICAL_VIRAL_TARGET** | <5% | Any | Essential for virus (may affect host) |
| **NON_SELECTIVE_TOXIC** | <50% | <80% | Would harm host |
| **NOT_A_TARGET** | >50% | Any | Minimal effect on virus |

**Key Results from Our Analysis**:
- **2 Selective Gene Targets**: PGM3 (selectivity: 100.0), GNPNAT1 (selectivity: 97.4)
- **7 Selective Reaction Targets**: UAGDP, ACGAMPM, GK1, ACGAM6PS, GF6PTA, METAT, AHCi
- **Target Pathway**: Hexosamine biosynthesis (UDP-GlcNAc synthesis)

### Core Modules (`src/`)

#### `config.py`
**Purpose**: Central configuration file with all constants and parameters.

**Key Contents**:
- File paths (genome, protein, model directories)
- HMPV protein copy numbers (from literature)
- Virion parameters (diameter, morphology)
- Lipid composition
- VBOF reaction parameters
- Metabolite ID mappings (BiGG IDs)

**Usage**: Import constants from this module:
```python
from src.config import HMPV_COPY_NUMBERS, VIRION_DIAMETER_NM, VBOF_REACTION_ID
```

---

#### `genome_analyzer.py`
**Purpose**: Analyzes HMPV genome sequence and calculates nucleotide stoichiometry.

**Key Functions**:
- `load_genome(path)`: Loads genome from FASTA file
- `parse_gff_annotations(path)`: Parses GFF annotation file
- `calculate_genome_stoichiometry(genome, copies_per_virion=1)`: Calculates nucleotide requirements

**Output**: Dictionary with BiGG metabolite IDs and coefficients:
```python
{
    'atp_c': -5103,  # ATP consumed
    'utp_c': -3402,  # UTP consumed
    'gtp_c': -2538,  # GTP consumed
    'ctp_c': -2307,  # CTP consumed
    'ppi_c': 13350   # Pyrophosphate produced
}
```

---

#### `protein_analyzer.py`
**Purpose**: Analyzes viral protein sequences and calculates amino acid stoichiometry.

**Key Functions**:
- `load_proteins(path)`: Loads proteins from FASTA file
- `calculate_protein_stoichiometry(proteins, copy_numbers)`: Calculates total amino acid requirements
- `set_copy_numbers(proteins, copy_numbers)`: Sets copy numbers for proteins

**Output**: Dictionary with amino acid BiGG IDs and coefficients:
```python
{
    'ala__L_c': -141190,
    'arg__L_c': -60070,
    # ... all 20 amino acids
}
```

---

#### `glycan_analyzer.py`
**Purpose**: Calculates glycan requirements for glycosylated proteins (F and G).

**Key Functions**:
- `calculate_f_protein_glycans(copy_number)`: Calculates F protein N-glycans
- `calculate_g_protein_glycans(copy_number)`: Calculates G protein N- and O-glycans
- `calculate_glycan_stoichiometry(f_copy_number, g_copy_number)`: Total glycan requirements

**Sources**:
- F protein: Viswanathan et al. (2011) - 3 N-linked sites
- G protein: Thammawat et al. (2008) - 5 N-linked + ~45 O-linked sites

**Output**: Dictionary with nucleotide sugar donors:
```python
{
    'uacgam_c': -7400,      # UDP-N-acetylglucosamine
    'gdpmann_c': -5550,     # GDP-mannose
    'udpgal_c': -11012,     # UDP-galactose
    # ... more glycan precursors
}
```

---

#### `vbof_builder.py`
**Purpose**: Assembles all components into complete VBOF.

**Key Functions**:
- `calculate_energy_requirements(genome_length, total_amino_acids, num_proteins)`: ATP/GTP costs
- `calculate_lipid_stoichiometry(virion_diameter_nm)`: Envelope lipid requirements
- `build_vbof(...)`: Main function to build complete VBOF

**Energy Costs**:
- RNA polymerization: ~2 ATP per nucleotide
- Translation: ~4 ATP per amino acid (activation + elongation)

**Lipid Calculation**:
- Calculates virion surface area from diameter
- Estimates lipid molecules based on packing density (0.65 nm² per lipid)
- Uses paramyxovirus lipid composition

**Output**: `VBOFComponents` dataclass with:
- `genome_stoichiometry`: Nucleotide requirements
- `protein_stoichiometry`: Amino acid requirements
- `energy_stoichiometry`: ATP/GTP requirements
- `lipid_stoichiometry`: Lipid requirements
- `glycan_stoichiometry`: Glycan requirements
- `combined_stoichiometry`: All combined
- `metadata`: Summary information

---

#### `model_integration.py`
**Purpose**: Integrates VBOF into host metabolic models using COBRApy.

**Key Functions**:
- `load_host_model(path)`: Loads host model from SBML file
- `map_vbof_metabolites(model, vbof_stoichiometry)`: Maps VBOF IDs to model IDs
- `create_vbof_reaction(model, vbof_stoichiometry)`: Creates COBRApy reaction
- `integrate_vbof(model, vbof_stoichiometry)`: Main integration function
- `validate_integrated_model(model)`: Validates integrated model

**Metabolite Mapping**:
- Attempts exact ID match first
- Suggests similar metabolites if not found
- Skips unmapped metabolites (with warning)

**VBOF Reaction**:
- Reaction ID: `HMPV_VBOF`
- Lower bound: 0 (irreversible)
- Upper bound: 1000
- Subsystem: "Viral Replication"

---

#### `exceptions.py`
**Purpose**: Custom exception classes for error handling.

**Exception Classes**:
- `HMPVModelError`: Base exception class
- `MissingGenomeError`: Genome file not found
- `MissingProteinDataError`: Protein data missing
- `MissingCopyNumberError`: Copy number not specified
- `EnergyDataMissingError`: Energy calculation failed
- `IDMappingError`: Metabolite ID mapping failed
- `ModelIntegrationError`: Model integration failed
- `VBOFConstructionError`: VBOF building failed

---

### Utility Scripts

#### `clean_host_model.py`
**Purpose**: Cleans host metabolic model (removes viral reactions, fixes issues).

**Usage**:
```bash
python clean_host_model.py
```

---

#### `generate_report.py`
**Purpose**: Generates comprehensive PDF report from analysis results.

**Usage**:
```bash
python generate_report.py
```

---

## Output Files

### VBOF Files (in `output_HBEC/`)

- **`hmpv_vbof.json`**: Raw VBOF with molecule counts
- **`hmpv_vbof_normalized.json`**: Normalized VBOF for FBA
- **`hmpv_vbof_summary.txt`**: Human-readable summary

### Integrated Model

- **`_with_HMPV_VBOF.xml`**: SBML file with integrated VBOF

### Dual-Objective Analysis (in `output_HBEC/dual_objective_analysis/`)

**Gene-Level Results:**
- **`host_growth_knockout_results.csv`**: 1,390 gene knockouts (host objective)
- **`virus_growth_knockout_results.csv`**: 1,390 gene knockouts (virus objective)
- **`merged_knockout_comparison.csv`**: Combined comparison with selectivity scores
- **`selective_antiviral_targets.csv`**: **KEY RESULT** - 2 selective targets
- **`critical_viral_targets.csv`**: 29 critical viral targets

**Reaction-Level Results:**
- **`host_reaction_knockout_results.csv`**: 1,769 reaction knockouts (host objective)
- **`virus_reaction_knockout_results.csv`**: Reaction knockouts (virus objective)
- **`merged_reaction_knockout_comparison.csv`**: Combined reaction comparison
- **`selective_reaction_targets.csv`**: **KEY RESULT** - 7 selective reaction targets
- **`critical_reaction_targets.csv`**: 58 critical viral reactions
- **`reaction_subsystem_essentiality.csv`**: Pathway-level analysis

**Reports:**
- **`dual_objective_report.txt`**: Comprehensive analysis report

### Figures (in `output_HBEC/figures/`)

- **`figure1_genome_map.png/svg`**: HMPV genome organization
- **`figure2_vbof_composition.png/svg`**: VBOF metabolite composition
- **`figure3_protein_copy_numbers.png/svg`**: Protein copy numbers per virion
- **`figure4_dual_objective_scatter.png/svg`**: Host vs virus growth scatter plot
- **`figure5_antiviral_targets.png/svg`**: Selective target visualization
- **`figure7_venn_diagram.png/svg`**: Target overlap analysis
- **`figure8_workflow.png/svg`**: Computational pipeline

### Research Documents (in project root)

- **`RESEARCH_PAPER.md`**: Complete research manuscript
- **`RESEARCH_HIGHLIGHTS.md`**: Key findings summary

### Sensitivity Analysis (if run)

- **`scenario_summary.csv`**: Summary of all parameter combinations
- **`robust_gene_targets.csv`**: Genes ranked by robustness score
- **`robust_reaction_targets.csv`**: Reactions ranked by robustness score
- **`comparison_report.txt`**: Cross-scenario comparison

---

## Configuration

### Key Parameters in `src/config.py`

#### Protein Copy Numbers
```python
HMPV_COPY_NUMBERS = {
    'N': 1900,    # Nucleoprotein
    'P': 300,     # Phosphoprotein
    'L': 30,      # RNA polymerase
    'M': 2000,    # Matrix protein
    'F': 350,     # Fusion protein
    'G': 250,     # Attachment glycoprotein
    'SH': 50,     # Small hydrophobic protein
    'M2-1': 100,  # Matrix protein 2-1
    'M2-2': 30,   # Matrix protein 2-2
}
```

#### Virion Parameters
```python
VIRION_DIAMETER_NM = 200.0  # Representative diameter (nm)
VIRION_DIAMETER_RANGE = (150.0, 600.0)  # Full range (pleomorphic)
```

#### Lipid Composition
```python
LIPID_COMPOSITION = {
    'pc_hs': {'fraction': 0.45},  # Phosphatidylcholine
    'pe_hs': {'fraction': 0.25},  # Phosphatidylethanolamine
    'ps_hs': {'fraction': 0.05},  # Phosphatidylserine
    'sphmyln_hs': {'fraction': 0.15},  # Sphingomyelin
    'chsterol': {'fraction': 0.10},  # Cholesterol
}
```

---

## References

### Key Papers

1. **Aller et al. (2018)** - "Integrated human-virus metabolic stoichiometric modelling predicts host-based antiviral targets"
   - *J R Soc Interface* 15:20180442
   - DOI: 10.1098/rsif.2018.0442

2. **Renz et al. (2021)** - "Genome-Scale Metabolic Model of Infection with SARS-CoV-2 Mutants"
   - *Genes* 12(6):796
   - DOI: 10.3390/genes12060796

3. **Viswanathan et al. (2011)** - F protein glycosylation
   - *J Gen Virol* 92:1580-1584
   - DOI: 10.1099/vir.0.030049-0

4. **Thammawat et al. (2008)** - G protein glycosylation
   - *J Virol* 82:10022-10034
   - DOI: 10.1128/JVI.01287-06

### Software Tools

- **COBRApy**: Constraint-based metabolic modeling
  - https://cobrapy.readthedocs.io/
  
- **BiGG Models**: Metabolic model database
  - http://bigg.ucsd.edu/

### Data Sources

- **HMPV Genome**: NCBI GenBank (GCF_002815375.1)
- **Host Model**: iHsaEC21 (Human epithelial cell model)

---

## Troubleshooting

### Common Issues

1. **"Baseline flux is zero"**
   - **Cause**: VBOF not normalized or model infeasible
   - **Solution**: Ensure `normalize_vbof.py` was run before integration

2. **"Metabolite not found"**
   - **Cause**: ID mismatch between VBOF and host model
   - **Solution**: Check metabolite mappings in `model_integration.py`

3. **"GPR.__init__() error"**
   - **Cause**: COBRApy version compatibility issue
   - **Solution**: Use fresh model loading (already fixed in sensitivity analysis)

4. **"No scenarios completed"**
   - **Cause**: All scenarios failing during processing
   - **Solution**: Check logs for specific errors, verify model file exists

---

## Contributing

This project follows strict coding standards:

- **PEP8** compliance
- **Type hints** for all functions
- **Docstrings** for all functions and classes
- **Modular architecture** - no monolithic scripts
- **Error handling** with custom exceptions
- **Logging** instead of print statements

---

## License

Copyright (c) 2025 Syed Mushahid Hussain and MSBI Dynamics - Justus-Liebig-Universität Gießen . All rights reserved.
<!--
1. License Grant
   This source code and associated data are provided solely for **personal review or educational purposes**.
   Use of this code is **strictly limited to non-commercial, non-distributable purposes**. 
   No part of this work may be copied, modified, shared, or incorporated into other work
   without the **explicit written permission** of the copyright holder.

2. Restrictions
   a. You may not distribute, upload, or share this project in any public or private repository.  
   b. You may not use this project in any publication, research, or commercial activity without written consent.  
   c. You may not reverse engineer or attempt to derive any part of the project for unauthorized purposes.  

3. Ownership
   All rights, including intellectual property rights, remain with the author. This license
   does not grant any rights to use the work in ways other than those explicitly permitted above.

4. Termination
   Any violation of this license will automatically terminate your rights to access or use the project
   and may result in legal action.

5. Disclaimer
   The project is provided "as-is" for informational purposes only. The author is not liable
   for any damages or misuse of the work.

By accessing or using this project, you acknowledge and agree to comply with all terms of this license.
-->

---

## Contact

**Author**: Syed Mushahid Hussain  
**Project**: HMPV Metabolic Modeling

For questions or issues, please refer to the documentation or check the code comments or email at hussains@students.uni-marburg.de.

---

## Version History

- **v3.0.0** (February 2026):
  - **NEW**: Dual-objective knockout analysis framework
  - **NEW**: Selective antiviral target identification (PGM3, GNPNAT1)
  - **NEW**: Reaction-level knockout analysis
  - **NEW**: Comprehensive research paper and highlights
  - **NEW**: Figure generation pipeline
  - Integration with HBEC metabolic model
  - Improved target classification system

- **v2.0.0** (December 2025):
  - VBOF construction pipeline
  - Model integration
  - Single-objective antiviral target analysis
  - Sensitivity analysis

---

*Last updated: February 2026*

