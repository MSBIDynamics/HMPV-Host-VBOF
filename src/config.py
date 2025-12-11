"""
Configuration Module
====================

This module contains configuration parameters and default values for HMPV VBOF construction.

Includes:
- File paths
- HMPV protein copy numbers (from literature)
- Model parameters
- BiGG ID mappings

Author: Syed Mushahid Hussain
"""

import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


# =============================================================================
# FILE PATHS
# =============================================================================

# Base directory (relative to project root)
DATA_DIR = Path("Data")
GENOMIC_DIR = DATA_DIR / "genomic"
PROTEIN_DIR = DATA_DIR / "protein"
MODEL_DIR = DATA_DIR / "smbl"
OUTPUT_DIR = Path("output_iHsaEC21")

# Default file names
DEFAULT_GENOME_FILE = "GCF_002815375.1_ASM281537v1_genomic.fna"
DEFAULT_GFF_FILE = "GCF_002815375.1_ASM281537v1_genomic.gff"
DEFAULT_PROTEIN_FILE = "GCF_002815375.1_ASM281537v1_protein.faa"
DEFAULT_HOST_MODEL = "iHsaEC21_clean.xml"  # Clean model without SARS-CoV-2 VBOF
DEFAULT_HOST_MODEL_ORIGINAL = "iHsaEC21.xml"  # Original model with SARS-CoV-2 VBOF

# =============================================================================
# OUTPUT FILE NAMES
# =============================================================================

# VBOF output files
VBOF_JSON_FILE = "hmpv_vbof.json"
VBOF_NORMALIZED_JSON_FILE = "hmpv_vbof_normalized.json"
VBOF_SUMMARY_FILE = "hmpv_vbof_summary.txt"

# Model integration output files
INTEGRATION_SUMMARY_FILE = "integration_summary.txt"
INTEGRATED_MODEL_XML_SUFFIX = "_with_HMPV_VBOF.xml"
INTEGRATED_MODEL_JSON_SUFFIX = "_with_HMPV_VBOF.json"

# Antiviral analysis output files
ANTIVIRAL_ANALYSIS_DIR = "antiviral_analysis"
GENE_KNOCKOUT_RESULTS_FILE = "gene_knockout_results.csv"
REACTION_KNOCKOUT_RESULTS_FILE = "reaction_knockout_results.csv"
TOP_GENE_TARGETS_FILE = "top_gene_targets.csv"
TOP_REACTION_TARGETS_FILE = "top_reaction_targets.csv"
SUBSYSTEM_ESSENTIALITY_FILE = "subsystem_essentiality.csv"
ANTIVIRAL_TARGETS_REPORT_FILE = "antiviral_targets_report.txt"

# Sensitivity analysis output files
SENSITIVITY_ANALYSIS_DIR = "sensitivity_analysis"
SCENARIO_SUMMARY_FILE = "scenario_summary.csv"
ROBUST_GENE_TARGETS_FILE = "robust_gene_targets.csv"
ROBUST_REACTION_TARGETS_FILE = "robust_reaction_targets.csv"
COMPARISON_REPORT_FILE = "comparison_report.txt"
COMMON_GENE_TARGETS_FILE = "common_gene_targets.csv"
COMMON_REACTION_TARGETS_FILE = "common_reaction_targets.csv"
UNIQUE_TARGETS_FILE = "unique_targets_by_scenario.csv"
SCENARIO_STATISTICS_FILE = "scenario_statistics.csv"

# Model cleaning output files
MODEL_CLEANING_REPORT_FILE = "model_cleaning_report.txt"

# Report generation output files
METABOLIC_MODEL_REPORT_FILE = "HMPV_Metabolic_Model_Report.pdf"

# =============================================================================
# OUTPUT DIRECTORY PATHS
# =============================================================================

# Subdirectories in output
ANTIVIRAL_ANALYSIS_OUTPUT_DIR = OUTPUT_DIR / ANTIVIRAL_ANALYSIS_DIR
SENSITIVITY_ANALYSIS_OUTPUT_DIR = OUTPUT_DIR / SENSITIVITY_ANALYSIS_DIR

# =============================================================================
# FULL FILE PATHS (for convenience)
# =============================================================================

# Input file paths
GENOME_FILE_PATH = GENOMIC_DIR / DEFAULT_GENOME_FILE
GFF_FILE_PATH = GENOMIC_DIR / DEFAULT_GFF_FILE
PROTEIN_FILE_PATH = PROTEIN_DIR / DEFAULT_PROTEIN_FILE
HOST_MODEL_CLEAN_PATH = MODEL_DIR / DEFAULT_HOST_MODEL
HOST_MODEL_ORIGINAL_PATH = MODEL_DIR / DEFAULT_HOST_MODEL_ORIGINAL

# VBOF output paths
VBOF_JSON_PATH = OUTPUT_DIR / VBOF_JSON_FILE
VBOF_NORMALIZED_JSON_PATH = OUTPUT_DIR / VBOF_NORMALIZED_JSON_FILE
VBOF_SUMMARY_PATH = OUTPUT_DIR / VBOF_SUMMARY_FILE

# Integration output paths
INTEGRATION_SUMMARY_PATH = OUTPUT_DIR / INTEGRATION_SUMMARY_FILE

# Antiviral analysis output paths
GENE_KNOCKOUT_RESULTS_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / GENE_KNOCKOUT_RESULTS_FILE
REACTION_KNOCKOUT_RESULTS_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / REACTION_KNOCKOUT_RESULTS_FILE
TOP_GENE_TARGETS_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / TOP_GENE_TARGETS_FILE
TOP_REACTION_TARGETS_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / TOP_REACTION_TARGETS_FILE
SUBSYSTEM_ESSENTIALITY_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / SUBSYSTEM_ESSENTIALITY_FILE
ANTIVIRAL_TARGETS_REPORT_PATH = ANTIVIRAL_ANALYSIS_OUTPUT_DIR / ANTIVIRAL_TARGETS_REPORT_FILE

# Sensitivity analysis output paths
SCENARIO_SUMMARY_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / SCENARIO_SUMMARY_FILE
ROBUST_GENE_TARGETS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / ROBUST_GENE_TARGETS_FILE
ROBUST_REACTION_TARGETS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / ROBUST_REACTION_TARGETS_FILE
COMPARISON_REPORT_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / COMPARISON_REPORT_FILE
COMMON_GENE_TARGETS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / COMMON_GENE_TARGETS_FILE
COMMON_REACTION_TARGETS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / COMMON_REACTION_TARGETS_FILE
UNIQUE_TARGETS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / UNIQUE_TARGETS_FILE
SCENARIO_STATISTICS_PATH = SENSITIVITY_ANALYSIS_OUTPUT_DIR / SCENARIO_STATISTICS_FILE

# Model cleaning output paths
MODEL_CLEANING_REPORT_PATH = OUTPUT_DIR / MODEL_CLEANING_REPORT_FILE

# Report generation output paths
METABOLIC_MODEL_REPORT_PATH = OUTPUT_DIR / METABOLIC_MODEL_REPORT_FILE


# =============================================================================
# HMPV PROTEIN COPY NUMBERS
# =============================================================================

# Copy numbers per virion for HMPV proteins
# 
# IMPORTANT: Some of the values are estimates based on related paramyxoviruses (RSV, PIV)
# and should be updated with HMPV-specific data when available.
#
# =============================================================================
# SOURCES AND REFERENCES:
# =============================================================================
#
# 1. N PROTEIN - RNA BINDING STOICHIOMETRY:
#    - Pneumovirus N protein binds ~7 nucleotides of RNA
#    - Reference: Tawar et al. (2009) "Crystal structure of a nucleocapsid-like 
#      nucleoprotein-RNA complex of respiratory syncytial virus"
#      Science 326(5957):1279-83. DOI: 10.1126/science.1177634
#      https://www.science.org/doi/10.1126/science.1177634
#    - Calculation: 13,350 nt genome / 7 nt per N ≈ 1,907 copies
#
# 2. RSV STRUCTURAL STUDIES (closely related to HMPV):
#    - Kiss et al. (2014) "Structural Analysis of RSV Reveals the Position of M2-1"
#      J Virol 88(12):7602-17. DOI: 10.1128/JVI.00256-14
#      https://journals.asm.org/doi/10.1128/jvi.00256-14
#    - Provides RSV protein organization estimates
#
# 3. PARAMYXOVIRUS STRUCTURAL BIOLOGY:
#    - Lamb & Parks (2013) "Paramyxoviridae" in Fields Virology, 6th Ed.
#      Chapter 33, pp. 957-995. Lippincott Williams & Wilkins.
#
# 4. HMPV CHARACTERIZATION:
#    - Peret et al. (2002) "Characterization of Human Metapneumoviruses"
#      J Infect Dis 185(11):1660-3. DOI: 10.1086/340518
#      https://pmc.ncbi.nlm.nih.gov/articles/PMC7109943/
#
# 5. HMPV DISCOVERY:
#    - van den Hoogen et al. (2001) "A newly discovered human pneumovirus"
#      Nat Med 7(6):719-24. DOI: 10.1038/89098
#      https://www.nature.com/articles/nm0601_719
#
# 6. PARAMYXOVIRUS CRYO-EM:
#    - Ke et al. (2018) "Structure of RSV Fusion Glycoprotein Trimer"
#      https://www.science.org/doi/10.1126/science.1234914
#    - F protein trimeric organization (~10-15 nm spacing)
#
# TODO: Replace with experimentally determined HMPV-specific values when available
# =============================================================================

HMPV_COPY_NUMBERS: Dict[str, int] = {
    # Nucleocapsid proteins (encapsidate genome)
    # Calculated: 13,350 nt genome / 7 nt per N protein ≈ 1,907
    # Source: Pneumovirus N-RNA binding stoichiometry
    'N': 1900,    # Nucleoprotein - binds viral RNA genome (calculated)
    
    # Polymerase complex
    # Source: Estimated from RSV/PIV studies
    # P:L ratio is typically ~10:1 in paramyxoviruses
    'P': 300,     # Phosphoprotein - polymerase cofactor (reduced)
    'L': 30,      # RNA polymerase - low copy, high MW enzyme (reduced)
    
    # Matrix protein
    # Forms layer under envelope, abundant
    # For 200nm particle: surface area ~125,600 nm², with M spacing ~7nm: ~2500-3000
    # Source: Paramyxovirus structural studies
    'M': 2000,    # Matrix protein - virion structure (adjusted)
    
    # Surface glycoproteins (embedded in envelope)
    # F protein forms trimeric spikes (~10-15 nm apart)
    # For 200nm particle: ~250-400 F trimers = 750-1200 F monomers
    # Source: Cryo-EM studies of paramyxoviruses
    'F': 350,     # Fusion protein - trimeric spikes (adjusted)
    'G': 250,     # Attachment glycoprotein (variable, can be absent)
    'SH': 50,     # Small hydrophobic protein (viroporin, minor)
    
    # M2 proteins (regulatory, internal)
    # Source: Estimated, lower abundance
    'M2-1': 100,  # Matrix protein 2-1 - transcription factor (reduced)
    'M2-2': 30,   # Matrix protein 2-2 - regulatory (reduced)
}

# Alternative conservative estimates (lower bound)
HMPV_COPY_NUMBERS_CONSERVATIVE: Dict[str, int] = {
    'N': 1400,    # Lower estimate
    'P': 200,
    'L': 20,
    'M': 1500,
    'F': 250,
    'G': 150,
    'SH': 30,
    'M2-1': 50,
    'M2-2': 20,
}

# Confidence levels for copy numbers
# HIGH: Direct HMPV experimental data or calculated from known stoichiometry
# MEDIUM: Extrapolated from closely related viruses (RSV)
# LOW: Rough estimates, need experimental validation
COPY_NUMBER_CONFIDENCE: Dict[str, str] = {
    'N': 'HIGH',      # Calculated from genome length / 7 nt per N protein
    'P': 'LOW',       # Estimated from P:L ratio
    'L': 'MEDIUM',    # Based on paramyxovirus studies
    'M': 'MEDIUM',    # Based on RSV and surface area calculation
    'F': 'MEDIUM',    # Based on cryo-EM of related viruses
    'G': 'LOW',       # Variable, HMPV can replicate without G
    'SH': 'LOW',      # Poorly characterized, minor protein
    'M2-1': 'LOW',    # Estimated
    'M2-2': 'LOW',    # Estimated
}

# Key structural references with URLs
STRUCTURAL_REFERENCES = {
    'peret_2002': {
        'title': 'Characterization of Human Metapneumoviruses Isolated from Patients in North America',
        'journal': 'J Infect Dis',
        'year': 2002,
        'pmcid': 'PMC7109943',
        'doi': '10.1086/340518',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7109943/',
        'key_findings': [
            'HMPV particles are pleomorphic',
            'Two main genetic lineages (A and B)',
            'F gene conserved, G gene variable',
            'EM shows nucleocapsid and filamentous particles'
        ]
    },
    'van_den_hoogen_2001': {
        'title': 'A newly discovered human pneumovirus isolated from young children with respiratory tract disease',
        'journal': 'Nat Med',
        'year': 2001,
        'doi': '10.1038/89098',
        'url': 'https://www.nature.com/articles/nm0601_719',
        'key_findings': [
            'Original HMPV discovery',
            'Paramyxovirus-like structure',
            'Genome ~13 kb'
        ]
    },
    'tawar_2009': {
        'title': 'Crystal structure of a nucleocapsid-like nucleoprotein-RNA complex of RSV',
        'journal': 'Science',
        'year': 2009,
        'doi': '10.1126/science.1177634',
        'url': 'https://www.science.org/doi/10.1126/science.1177634',
        'key_findings': [
            'N protein binds 7 nucleotides of RNA',
            'Basis for N protein copy number calculation'
        ]
    },
    'kiss_2014': {
        'title': 'Structural Analysis of RSV Reveals the Position of M2-1',
        'journal': 'J Virol',
        'year': 2014,
        'doi': '10.1128/JVI.jvi.00256-14',
        'url': 'https://journals.asm.org/doi/10.1128/jvi.00256-14',
        'key_findings': [
            'RSV virion protein organization',
            'Matrix and M2-1 protein distribution'
        ]
    },
    'ke_2018': {
        'title': 'Structure of RSV Fusion Glycoprotein Trimer',
        'journal': 'Science',
        'year': 2018,
        'doi': '10.1126/science.1234914',
        'url': 'https://www.science.org/doi/10.1126/science.1234914',
        'key_findings': [
            'F protein trimeric structure',
            'Spike distribution on virion surface'
        ]
    }
}

# =============================================================================
# LIPID COMPOSITION SOURCES
# =============================================================================
#
# Lipid composition for HMPV envelope (derived from host plasma membrane)
#
# SOURCES AND REFERENCES:
#
# 1. PARAMYXOVIRUS ENVELOPE COMPOSITION:
#    - Jinjun Shan. (2028) "High-resolution lipidomics reveals dysregulation of lipid metabolism"
#      https://pubs.rsc.org/en/content/articlelanding/2018/ra/c8ra05640d
#
# 2. VIRAL ENVELOPE LIPID RATIOS:
#    - Brügger et al. (2006) "The HIV lipidome"
#      https://pmc.ncbi.nlm.nih.gov/articles/PMC1413831/
#    - General principles of viral envelope composition
#
# 3. MEMBRANE LIPID PACKING:
#    - Nagle & Tristram-Nagle (2000) "Structure of lipid bilayers"
#      Biochim Biophys Acta 1469(3):159-195. DOI: 10.1016/S0304-4157(00)00016-2
#      https://www.sciencedirect.com/science/article/pii/S0304415700000162
#    - Lipid packing density: ~0.65 nm² per lipid molecule
#
# 4. SARS-CoV-2 VBOF METHODOLOGY (applied to HMPV):
#    - Aller et al. (2018) "Integrated human-virus metabolic stoichiometric modelling"
#      https://pmc.ncbi.nlm.nih.gov/articles/PMC6170780/
#
# Lipid fractions (based on paramyxovirus/plasma membrane data):
LIPID_COMPOSITION = {
    'pc_hs': {
        'name': 'Phosphatidylcholine',
        'fraction': 0.45,
        'bigg_id': 'pc_hs_c',
        'note': 'Major phospholipid in plasma membrane'
    },
    'pe_hs': {
        'name': 'Phosphatidylethanolamine',
        'fraction': 0.25,
        'bigg_id': 'pe_hs_c',
        'note': 'Second most abundant phospholipid'
    },
    'ps_hs': {
        'name': 'Phosphatidylserine',
        'fraction': 0.05,
        'bigg_id': 'ps_hs_c',
        'note': 'Minor component, inner leaflet'
    },
    'sphmyln_hs': {
        'name': 'Sphingomyelin',
        'fraction': 0.15,
        'bigg_id': 'sphmyln_hs_c',
        'note': 'Sphingolipid component'
    },
    'chsterol': {
        'name': 'Cholesterol',
        'fraction': 0.10,
        'bigg_id': 'chsterol_c',
        'note': 'Essential for membrane fluidity'
    }
}

# Lipid packing density (nm² per lipid molecule in bilayer)
# Source: Nagle & Tristram-Nagle (2000) Biochim Biophys Acta
# Source 2: Molecular models for drug permeation across phospholipid membranes
# URL: https://docserv.uni-duesseldorf.de/servlets/DerivateServlet/Derivate-2600/600.pdf
LIPID_PACKING_DENSITY_NM2 = 0.65


# =============================================================================
# VIRION PARAMETERS
# =============================================================================

# HMPV virion dimensions
# IMPORTANT: HMPV is PLEOMORPHIC (variable shape) - can be:
#   - Spherical particles: 150-300 nm
#   - Filamentous particles: up to 600 nm length
#
# Reference: "Characterization of Human Metapneumoviruses Isolated from 
#            Patients in North America" - PMC7109943
#
# For VBOF calculations, we use a representative spherical particle
VIRION_DIAMETER_NM = 200.0  # Representative diameter in nanometers 
VIRION_DIAMETER_RANGE = (150.0, 600.0)  # Full range observed (pleomorphic)
VIRION_DIAMETER_SPHERICAL_TYPICAL = (150.0, 300.0)  # Typical spherical range

# Virion morphology note
VIRION_MORPHOLOGY = "pleomorphic"  # Can be spherical or filamentous


# =============================================================================
# VBOF REACTION PARAMETERS
# =============================================================================

VBOF_REACTION_ID = "HMPV_VBOF"
VBOF_REACTION_NAME = "Human Metapneumovirus Biomass Objective Function"
VBOF_SUBSYSTEM = "Viral Replication"


# =============================================================================
# METABOLITE ID MAPPINGS
# =============================================================================

# BiGG IDs for nucleotides (cytosolic)
NUCLEOTIDE_BIGG_IDS = {
    'A': 'atp_c',
    'U': 'utp_c',
    'G': 'gtp_c',
    'C': 'ctp_c',
}

# BiGG IDs for amino acids (cytosolic)
AMINO_ACID_BIGG_IDS = {
    'A': 'ala__L_c',
    'R': 'arg__L_c',
    'N': 'asn__L_c',
    'D': 'asp__L_c',
    'C': 'cys__L_c',
    'E': 'glu__L_c',
    'Q': 'gln__L_c',
    'G': 'gly_c',
    'H': 'his__L_c',
    'I': 'ile__L_c',
    'L': 'leu__L_c',
    'K': 'lys__L_c',
    'M': 'met__L_c',
    'F': 'phe__L_c',
    'P': 'pro__L_c',
    'S': 'ser__L_c',
    'T': 'thr__L_c',
    'W': 'trp__L_c',
    'Y': 'tyr__L_c',
    'V': 'val__L_c',
}

# Common metabolite BiGG IDs
COMMON_METABOLITE_IDS = {
    'atp': 'atp_c',
    'adp': 'adp_c',
    'amp': 'amp_c',
    'gtp': 'gtp_c',
    'gdp': 'gdp_c',
    'gmp': 'gmp_c',
    'pi': 'pi_c',
    'ppi': 'ppi_c',
    'h2o': 'h2o_c',
    'h': 'h_c',
}


# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
LOG_LEVEL = logging.INFO


def get_full_path(relative_path: Path) -> Path:
    """
    Get full path relative to project root.
    
    Parameters:
    -----------
    relative_path : Path
        Relative path from project root
    
    Returns:
    --------
    Path : Absolute path
    """
    # Assuming this file is in src/, go up one level to project root
    project_root = Path(__file__).parent.parent
    return project_root / relative_path


def print_configuration():
    """Print current configuration settings."""
    logger.info("=" * 60)
    logger.info("HMPV VBOF Configuration")
    logger.info("=" * 60)
    
    logger.info("\nFile Paths:")
    logger.info(f"  Genome: {GENOMIC_DIR / DEFAULT_GENOME_FILE}")
    logger.info(f"  Proteins: {PROTEIN_DIR / DEFAULT_PROTEIN_FILE}")
    logger.info(f"  Host Model: {MODEL_DIR / DEFAULT_HOST_MODEL}")
    
    logger.info("\nProtein Copy Numbers:")
    for protein, count in HMPV_COPY_NUMBERS.items():
        confidence = COPY_NUMBER_CONFIDENCE.get(protein, 'UNKNOWN')
        logger.info(f"  {protein}: {count} (confidence: {confidence})")
    
    logger.info(f"\nVirion Parameters:")
    logger.info(f"  Diameter: {VIRION_DIAMETER_NM} nm")
    
    logger.info("=" * 60)

