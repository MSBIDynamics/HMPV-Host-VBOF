"""
Configuration Module
====================

This module contains configuration parameters and default values for HMPV VBOF construction.

Includes:
- File paths
- Default HMPV protein copy numbers
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
MODEL_DIR = DATA_DIR / "sbml"
OUTPUT_DIR = Path("output_HBEC")  # Output directory for HBEC model

# Default file names
DEFAULT_GENOME_FILE = "GCF_002815375.1_ASM281537v1_genomic.fna"
DEFAULT_GFF_FILE = "GCF_002815375.1_ASM281537v1_genomic.gff"
DEFAULT_PROTEIN_FILE = "GCF_002815375.1_ASM281537v1_protein.faa"
DEFAULT_HOST_MODEL = "iHBEC_Recon3D_or.xml"  # iHBEC model with R_biomass_hbec, without SARS-CoV-2 VBOF
DEFAULT_HOST_MODEL_WITH_SARS = "iHBEC_Recon3D_or.xml"  # Original model with both R_biomass_hbec and R_VBOF (SARS-CoV-2)
DEFAULT_HOST_MODEL_ORIGINAL = "iHsaEC21.xml"  # Original iHsaEC21 model (alternative)


# =============================================================================
# HOST BOF REACTION PARAMETERS
# =============================================================================

HOST_BOF_REACTION_ID = "biomass_reaction"
HOST_BOF_REACTION_NAME = "Human Bronchial Epithelial Cell Biomass"

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

# Dual-objective analysis output files
DUAL_OBJECTIVE_ANALYSIS_DIR = "dual_objective_analysis"
HOST_GROWTH_RESULTS_FILE = "host_growth_knockout_results.csv"
VIRUS_GROWTH_RESULTS_FILE = "virus_growth_knockout_results.csv"
MERGED_KNOCKOUT_RESULTS_FILE = "merged_knockout_comparison.csv"
SELECTIVE_TARGETS_FILE = "selective_antiviral_targets.csv"
CRITICAL_VIRAL_TARGETS_FILE = "critical_viral_targets.csv"
COMBINED_OBJECTIVE_RESULTS_FILE = "combined_objective_results.csv"
DUAL_OBJECTIVE_REPORT_FILE = "dual_objective_report.txt"
# Reaction knockout files for dual-objective analysis
HOST_REACTION_KNOCKOUT_RESULTS_FILE = "host_reaction_knockout_results.csv"
VIRUS_REACTION_KNOCKOUT_RESULTS_FILE = "virus_reaction_knockout_results.csv"
MERGED_REACTION_KNOCKOUT_RESULTS_FILE = "merged_reaction_knockout_comparison.csv"
SELECTIVE_REACTION_TARGETS_FILE = "selective_reaction_targets.csv"
CRITICAL_REACTION_TARGETS_FILE = "critical_reaction_targets.csv"
REACTION_SUBSYSTEM_ESSENTIALITY_FILE = "reaction_subsystem_essentiality.csv"
SELECTIVE_TARGETS_FVA_FILE = "selective_targets_fva.csv"



# Multi-variant comparison output files
MULTI_VARIANT_COMPARISON_DIR = "multi_variant_comparison"
VARIANT_COMPARISON_SUMMARY_FILE = "variant_comparison_summary.csv"
VARIANT_COMPARISON_REPORT_FILE = "variant_comparison_report.txt"
MERGED_SELECTIVE_GENE_FILE = "merged_selective_gene_targets.csv"
MERGED_SELECTIVE_RXN_FILE = "merged_selective_rxn_targets.csv"
MERGED_FVA_FILE = "merged_fva_results.csv"


# =============================================================================
# OUTPUT DIRECTORY PATHS
# =============================================================================

# Subdirectories in output
DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR = OUTPUT_DIR / DUAL_OBJECTIVE_ANALYSIS_DIR

# =============================================================================
# FULL FILE PATHS (for convenience)
# =============================================================================

# Input file paths
GENOME_FILE_PATH = GENOMIC_DIR / DEFAULT_GENOME_FILE
GFF_FILE_PATH = GENOMIC_DIR / DEFAULT_GFF_FILE
PROTEIN_FILE_PATH = PROTEIN_DIR / DEFAULT_PROTEIN_FILE
HOST_MODEL_PATH = MODEL_DIR / DEFAULT_HOST_MODEL  # Primary model (model_clean.xml)
HOST_MODEL_CLEAN_PATH = MODEL_DIR / DEFAULT_HOST_MODEL  # Alias for backward compatibility
HOST_MODEL_WITH_SARS_PATH = MODEL_DIR / DEFAULT_HOST_MODEL_WITH_SARS  # Model with SARS-CoV-2 VBOF
HOST_MODEL_ORIGINAL_PATH = MODEL_DIR / DEFAULT_HOST_MODEL_ORIGINAL

# VBOF output paths
VBOF_JSON_PATH = OUTPUT_DIR / VBOF_JSON_FILE
VBOF_NORMALIZED_JSON_PATH = OUTPUT_DIR / VBOF_NORMALIZED_JSON_FILE
VBOF_SUMMARY_PATH = OUTPUT_DIR / VBOF_SUMMARY_FILE

# Integration output paths
INTEGRATION_SUMMARY_PATH = OUTPUT_DIR / INTEGRATION_SUMMARY_FILE

# Dual-objective analysis output paths
HOST_GROWTH_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / HOST_GROWTH_RESULTS_FILE
VIRUS_GROWTH_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / VIRUS_GROWTH_RESULTS_FILE
MERGED_KNOCKOUT_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / MERGED_KNOCKOUT_RESULTS_FILE
SELECTIVE_TARGETS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / SELECTIVE_TARGETS_FILE
CRITICAL_VIRAL_TARGETS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / CRITICAL_VIRAL_TARGETS_FILE
COMBINED_OBJECTIVE_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / COMBINED_OBJECTIVE_RESULTS_FILE
DUAL_OBJECTIVE_REPORT_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / DUAL_OBJECTIVE_REPORT_FILE
# Reaction knockout paths for dual-objective analysis
HOST_REACTION_KNOCKOUT_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / HOST_REACTION_KNOCKOUT_RESULTS_FILE
VIRUS_REACTION_KNOCKOUT_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / VIRUS_REACTION_KNOCKOUT_RESULTS_FILE
MERGED_REACTION_KNOCKOUT_RESULTS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / MERGED_REACTION_KNOCKOUT_RESULTS_FILE
SELECTIVE_REACTION_TARGETS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / SELECTIVE_REACTION_TARGETS_FILE
CRITICAL_REACTION_TARGETS_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / CRITICAL_REACTION_TARGETS_FILE
REACTION_SUBSYSTEM_ESSENTIALITY_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / REACTION_SUBSYSTEM_ESSENTIALITY_FILE
SELECTIVE_TARGETS_FVA_PATH = DUAL_OBJECTIVE_ANALYSIS_OUTPUT_DIR / SELECTIVE_TARGETS_FVA_FILE



# =============================================================================
# HMPV PROTEIN COPY NUMBERS
# =============================================================================

USE_CALCULATED_COPY_NUMBERS: bool = True


N_NUCLEOTIDES_PER_PROTOMER: int = 7    # nucleotides bound per N protomer
GENOME_COPIES_PER_VIRION: int = 1      # 1 for spherical particles; 2 for filamentous RSV

# --- M protein (geometric surface-area model) ---
# Source: Kiss et al. (2014)  (RSV cryo-ET, 24% M coverage)
#       Leyrat et al. (2014) PDB 4LP7 (HMPV M dimer convex-hull footprint)
M_MEMBRANE_COVERAGE_FRACTION: float = 0.24   # fraction of inner membrane covered by M layer
M_DIMER_FOOTPRINT_NM2: float = 35.8          # membrane-facing convex-hull area of M dimer (YZ plane, PDB 4LP7)

# --- Glycoproteins F, G, SH (hexagonal spike-packing model) ---
# Source: Walsh et al. (2015) — RSV spike width/spacing by EM
#         Conley et al. (2022) — RSV cryo-ET hexagonal arrays
#         McLellan (2013), Battles (2017) — F trimer; Thammawat (2008) — G monomer
#         Gan et al. (2012) — RSV SH pentamer (applied by analogy to HMPV SH)
SPIKE_WIDTH_NM: float = 11.5          # average glycoprotein spike width in nm
SPIKE_SPACING_MIN_NM: float = 6.0     # minimum inter-spike gap in nm
SPIKE_SPACING_AVG_NM: float = 8.0     # average inter-spike gap in nm
SPIKE_SPACING_MAX_NM: float = 10.0    # maximum inter-spike gap in nm
SPIKE_SPACING_MODE: str = 'average'   # which spacing to use: 'min', 'average', or 'max'

# F:G:SH stoichiometric ratio (based on transcription gradient + functional essentiality)
# Source: Schildgen et al. (2011) — transcription gradient; Biacchesi et al. (2004)
GLYCOPROTEIN_RATIO: Dict[str, int] = {'F': 5, 'G': 3, 'SH': 1}

# Oligomeric states per spike (subunits per spike position)
F_OLIGOMERIC_STATE: int = 3   # homotrimer  (McLellan 2013, Battles 2017)
G_OLIGOMERIC_STATE: int = 1   # monomer     (Thammawat 2008, Leyrat 2014)
SH_OLIGOMERIC_STATE: int = 5  # pentamer    (Gan 2012 — RSV SH analogy; Masante 2014)


HMPV_COPY_NUMBERS: Dict[str, int] = {
    # Nucleocapsid proteins (encapsidate genome)
    # Calculated: 13,350 nt genome / 7 nt per N protein ≈ 1,907
    # Source: Pneumovirus N-RNA binding stoichiometry
    'N': 1900,    # Nucleoprotein - binds viral RNA genome (calculated)
   
    'P': 300,     # From Sendai virus https://pdf.sciencedirectassets.com/272412/1-s2.0-S0042682200X01456/1-s2.0-S0042682296903591/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEBQaCXVzLWVhc3QtMSJIMEYCIQDeC0gSA%2BKrEkQOz%2Fe5i5O1zhX%2BmWOYJOpN12Q8H8TSLwIhAIak09xe6DgQF1b%2BuWkD7RaQzm4nZ0do29To21W7lRP1KrwFCN3%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBRoMMDU5MDAzNTQ2ODY1Igzdim2X5mUUiy2dhRQqkAV%2Bt8KxqrK0cEG7hjD4b9S832%2F5lmjOfM2rA0BIGoezzFN5fuBREV9Oak%2BXweefm3YnHN08NzdCN%2BnjTr3ok2WcIW1uKxefthb%2FTAkPBQlzP%2FI0oLcBB6gyrQti5vNjHKB8yFeEl1Aax%2B3HdiQ6TCPbXAfbGOoJbaPHohNWdk4QL4UL%2BUKHxVte6scEDGesGXGWMi1NAZ8a5NbqVUX8rzAqHbS0MsKrevc%2FrV3qSYzNktwNS8yh7vjSai9CnIIo8sYHspaIm9uf81FrOkwKXOmC4VXGY2ZsgD3kJLalvzsR50Yi4UO%2BmauK14Uh7DwhYexSCiQD9wQT%2FoM1Ir7f3I3egVPzNu5pdNsOuVhgpNe8U0ZKpmU5%2BFTAVaKTId0G22IQeQwlKxP43m6MGW9yDK9ulvNmXI%2BtPco%2FL6m7o4gpf%2Fs5ZVrENaAXeJgVpIqb3CqpTZfmlM45zZKjqzCtwBtjOTCIAeeyUlX9NsKH6XztP85uxu%2BDtTaxa35mHX%2FewjsacDIDyQiYyFc4%2BRk2A8Lbhe2lTIY1Q1%2F5ZqP1pFh92Z70DT3FbGReHdcx%2Bmt5ktyDffYZAJYHaO5isnH%2FodYp%2BXqBQs0MYMVlgYqfLyWGq7OwydD044588NDvqdNq%2FXTEvHupsuqIh9OYAeQRjRt5ialB6Jfh%2BLmwMoxaoE8j%2FL%2BjaPQph0WuRw9Ss9yEQG8gTBwf5lC0NgV4VMUz8JYzVmNw8P1QuPWyuiyJ461pi5sLwvmYg9arOISdASeF0H1KaqvfzUeSNCmUN1Vx5hk7Q1IxQdumoMIoNn09JlhsVjn5GHEgdnMUJ1i0QcDgJ8ToeVoVS0RA4Xr2JTcDAXILGQu8sb0x9j0OIjXWU9ZoCjDb9tHOBjqwAYTbJnhfYeKmdbpd%2Biy7lgYJ2kV1adTUangv1EmPBhpGxt%2B434z3199lYF7SJSpFgPcOVd13DLxdT1ernanWhduGaxNa7Lp5UW%2FPwQ0yAq6HIZjQAfwW6pcaqMo%2BcRa1%2FmEX4AZ8Q1tVgXaHLflG0HHuUHfQhFWiOC23XtkmoJDOjG%2FdkgSfKENDkuFLKjIpNOlLX8jjh3R4cUcKjlrITRfgRUE%2BRzmCRgD0%2FKGsawVa&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20260407T043126Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY4UGMCXQD%2F20260407%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=c2c0c9b59f16ef8b99e3205f30b4bb3114277a830d5792c6688fd973c334c079&hash=3162b9416debf0885161b80ce2594b9e38b16532a925dc56873a48762f376c79&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0042682296903591&tid=spdf-ecce83ef-2907-4a0c-b173-b1d973d19482&sid=3ad2fd778f23b3454f5bee5514227e4ec730gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=04015f0a5a57560150&rr=9e864d7a29a2290f&cc=de
    'L': 50,      # From Sendai virus
    
   
    'M': 180,    # Matrix protein - virion structure
    
    # Surface glycoproteins (embedded in envelope)
   
    'F': 693,     # Fusion protein - trimeric spikes 
    'G': 139,     # Attachment glycoprotein - monomeric spikes
    'SH': 230,    
    
    # M2 proteins (regulatory, internal)
    # Source: Estimated, lower abundance
    'M2-1': 200, 
    'M2-2': 50,  
}


# Lipid packing density (nm² per lipid molecule in bilayer)
# Source: Nagle & Tristram-Nagle (2000) Biochim Biophys Acta
# Source 2: Molecular models for drug permeation across phospholipid membranes
# URL: https://docserv.uni-duesseldorf.de/servlets/DerivateServlet/Derivate-2600/600.pdf
LIPID_PACKING_DENSITY_NM2: float = 0.65

# Lipid composition fractions (must sum to ≤ 1.0)
# Source: Barnes et al. (1987) — Sendai virus 
LIPID_FRACTIONS: Dict[str, float] = {
    'pchol_hs_c':   0.1865,   # Phosphatidylcholine       (18.65%)
    'pe_hs_c':      0.1340,   # Phosphatidylethanolamine  (13.40%)
    'ps_hs_c':      0.0600,   # Phosphatidylserine        ( 6.00%)
    'sphmyln_hs_c': 0.0440,   # Sphingomyelin             ( 4.40%)
    'clpn_hs_c':    0.0450,   # Cardiolipin               ( 4.50%)
    'pail_hs_c':    0.0275,   # Phosphatidylinositol      ( 2.75%)
    'lpe_hs_c':     0.0080,   # Lysophosphatidylethanolamine ( 0.80%)
    'chsterol_c':   0.9000,   # Cholesterol               (50.00%)
}


# =============================================================================
# VIRION PARAMETERS
# =============================================================================

VIRION_DIAMETER_NM = 209.0  # Default representative diameter (nm)
VIRION_DIAMETER_RANGE = (150.0, 600.0)  # Full range observed (pleomorphic)
VIRION_DIAMETER_SPHERICAL_TYPICAL = (150.0, 300.0)  # Typical spherical range

# Virion morphology note
VIRION_MORPHOLOGY = "pleomorphic"  # Can be spherical or filamentous

FILAMENTOUS_DIAMETER_NM: float = 120.0   # tube cross-section diameter (nm)
FILAMENTOUS_LENGTH_NM:   float = 600.0   # total tip-to-tip length (nm)


# =============================================================================
# GLYCAN PARAMETERS
# =============================================================================
#

F_PROTEIN_N_LINKED_SITES: int = 3         # sites N57, N172, N353 — all utilized
F_PROTEIN_SITE_UTILIZATION: float = 1.0   # 100 % occupancy confirmed

G_PROTEIN_N_LINKED_SITES: int = 4         # average across genotypes (range 2–5)
G_PROTEIN_N_LINKED_UTILIZATION: float = 0.8   # ~80 % occupancy
G_PROTEIN_O_LINKED_SITES: int = 26        # conservative estimate of occupied sites
G_PROTEIN_O_LINKED_UTILIZATION: float = 1.0  # ~50 % of potential Ser/Thr sites

N_GLYCAN_COMPLEX_COMPOSITION: Dict[str, int] = {
    'GlcNAc':  4,   # 2 core + 2 antennae
    'Man':     3,   # core trimannose
    'Gal':     2,   # terminal galactose on each antenna
    'Neu5Ac':  2,   # sialic acid caps
    'Fuc':     1,   # core fucose
}
N_LINKED_GLYCAN_HIGH_MANNOSE = {
        'GlcNAc': 2,    # Core chitobiose
        'Man': 9,       # High-mannose
}
# O-linked glycan monosaccharide composition — Core 1 sialyl-T antigen (mucin-type)

O_GLYCAN_CORE1_COMPOSITION: Dict[str, int] = {
    'GalNAc':  1,   # N-acetylgalactosamine (linked to Ser/Thr)
    'Gal':     1,   # galactose
    'Neu5Ac':  1,   # sialic acid cap
}

# ATP cost per glycosidic bond formation
GLYCOSYLATION_ATP_PER_N_GLYCAN: int = 2   # per N-glycan (OST reaction)
GLYCOSYLATION_ATP_PER_O_GLYCAN: int = 1   # per O-glycan (GalNAc-T reaction)


# =============================================================================
# VBOF REACTION PARAMETERS
# =============================================================================

VBOF_REACTION_ID = "HMPV_VBOF"
VBOF_REACTION_NAME = "Human Metapneumovirus Biomass Objective Function"
VBOF_SUBSYSTEM = "Viral Replication"



# =============================================================================
# VIRION DIAMETER VARIANTS  (for multi-variant pipeline)
# =============================================================================
#

VIRION_DIAMETER_VARIANTS: Dict[str, Dict] = {
    "spherical_small": {
        "morphology":    "spherical",
        "diameter_nm":   150.0,
        "length_nm":     None,
        "genome_copies": 1,
    },
    "spherical_typical": {
        "morphology":    "spherical",
        "diameter_nm":   209.0,
        "length_nm":     None,
        "genome_copies": 1,
    },
    "spherical_large": {
        "morphology":    "spherical",
        "diameter_nm":   300.0,
        "length_nm":     None,
        "genome_copies": 1,
    },
    "filamentous_small": {
        "morphology":    "filamentous",
        "diameter_nm":   45,   
        "length_nm":     200,      
        "genome_copies": 1,    
    },
    "filamentous_medium": {
        "morphology":    "filamentous",
        "diameter_nm":   62,  
        "length_nm":     282,      
        "genome_copies": 1,    
    },
    "filamentous_large": {
        "morphology":    "filamentous",
        "diameter_nm":   120,  
        "length_nm":     600,      
        "genome_copies": 1,    
    },
}


# =============================================================================
# DUAL-OBJECTIVE ANALYSIS PARAMETERS
# =============================================================================

# Default threshold configurations for target classification
# Values are expressed as fractions (0-1), where 1.0 = 100% of wild-type flux
DEFAULT_THRESHOLDS: Dict[str, Dict[str, float]] = {
    # Selective target: significantly hurts virus while sparing host
    'selective_target': {
        'virus_max': 0.5,   # Virus growth drops below 50%
        'host_min': 0.8     # Host growth stays above 80%
    },
    # Strict selective target: strong effect on virus, minimal effect on host
    'strict_selective': {
        'virus_max': 0.1,   # Virus growth drops below 10%
        'host_min': 0.9     # Host growth stays above 90%
    },
    # Lethal to virus: critical for virus survival
    'lethal_virus': {
        'virus_max': 0.05   # Virus growth drops below 5% (almost zero)
    }
}

# Classification thresholds for impact assessment
ESSENTIALITY_THRESHOLD = 0.01      # < 1% of max flux = lethal
SIGNIFICANT_THRESHOLD = 0.5         # < 50% of max flux = significant
MODERATE_THRESHOLD = 0.9            # < 90% of max flux = moderate


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

GLYCAN_PRECURSOR_BIGG_IDS: Dict[str, str] = {
    'GlcNAc':  'uacgam_c',    # UDP-N-acetyl-D-glucosamine
    'Man':     'gdpmann_c',   # GDP-mannose
    'Gal':     'udpgal_c',    # UDP-galactose
    'Fuc':     'gdpfuc_c',    # GDP-fucose (cytosol)
    'Neu5Ac':  'cmpacna_c',   # CMP-N-acetylneuraminate (sialic acid)
    'GalNAc':  'uacgam_c',  # UDP-N-acetyl-D-galactosamine
}

# Molecular weights (g/mol) for all VBOF metabolites
METABOLITE_MOLECULAR_WEIGHTS: Dict[str, float] = {
    # Genome nucleotides
    'atp_c': 507.18,
    'gtp_c': 523.18,
    'ctp_c': 483.16,
    'utp_c': 484.14,
    # Amino acids (all 20 standard)
    'ala__L_c': 89.09,
    'arg__L_c': 174.20,
    'asn__L_c': 132.12,
    'asp__L_c': 133.10,
    'cys__L_c': 121.16,
    'glu__L_c': 147.13,
    'gln__L_c': 146.15,
    'gly_c':    75.07,
    'his__L_c': 155.16,
    'ile__L_c': 131.17,
    'leu__L_c': 131.17,
    'lys__L_c': 146.19,
    'met__L_c': 149.21,
    'phe__L_c': 165.19,
    'pro__L_c': 115.13,
    'ser__L_c': 105.09,
    'thr__L_c': 119.12,
    'trp__L_c': 204.23,
    'tyr__L_c': 181.19,
    'val__L_c': 117.15,
    # Lipids 
    'pchol_hs_c':   760.0,    # Phosphatidylcholine
    'pe_hs_c':      720.0,    # Phosphatidylethanolamine
    'ps_hs_c':      780.0,    # Phosphatidylserine
    'sphmyln_hs_c': 730.0,    # Sphingomyelin
    'chsterol_c':   386.65,   # Cholesterol
    'clpn_hs_c':    1448.0,   # Cardiolipin
    'pail_hs_c':    878.0,    # Phosphatidylinositol
    'lpe_hs_c':     479.0,    # Lysophosphatidylethanolamine
    # Glycan precursors
    'uacgam_c':   607.35,   # UDP-N-acetylglucosamine
    'gdpmann_c':  605.34,   # GDP-mannose
    'udpgal_c':   566.30,   # UDP-galactose
    'cmpacna_c':  614.39,   # CMP-N-acetylneuraminic acid (sialic acid)
    'gdpfuc_c':   589.33,   # GDP-fucose
    'M_uacgam_c': 607.35,   # UDP-N-acetylgalactosamine
}

# Metabolites physically incorporated into the virion (used for mass calculations).
STRUCTURAL_METABOLITES: set = (
    set(AMINO_ACID_BIGG_IDS.values())           # all 20 amino acids
    | set(LIPID_FRACTIONS.keys())               # all lipid envelope species
    | set(GLYCAN_PRECURSOR_BIGG_IDS.values())   # all glycan precursors
    | {NUCLEOTIDE_BIGG_IDS['C'], NUCLEOTIDE_BIGG_IDS['U']}  # genome C and U only
)


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
        logger.info(f"  {protein}: {count}")
    
    logger.info(f"\nVirion Parameters:")
    logger.info(f"  Diameter: {VIRION_DIAMETER_NM} nm")
    
    logger.info("=" * 60)

