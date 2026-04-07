#!/usr/bin/env python3
"""
HMPV VBOF Construction Pipeline
===============================

This script builds the Viral Biomass Objective Function (VBOF) for Human Metapneumovirus (HMPV)
using the genome and protein data.

Usage:
------
    python build_vbof.py

Output:
-------
    - output/hmpv_vbof.json: Complete VBOF stoichiometry
    - output/hmpv_vbof_summary.txt: Human-readable summary
    - Console output with detailed logs

Author: Syed Mushahid Hussain
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime

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
    get_protein_summary,
    get_copy_numbers,
)
from src.vbof_builder import (
    build_vbof,
    get_vbof_summary,
    export_vbof_to_dict
)
from src.config import (
    GENOMIC_DIR,
    PROTEIN_DIR,
    OUTPUT_DIR,
    DEFAULT_GENOME_FILE,
    DEFAULT_GFF_FILE,
    DEFAULT_PROTEIN_FILE,
    VBOF_JSON_PATH,
    VBOF_SUMMARY_PATH,
    DEFAULT_GENOME_FILE,
    DEFAULT_GFF_FILE,
    DEFAULT_PROTEIN_FILE,
    HMPV_COPY_NUMBERS,
    COPY_NUMBER_CONFIDENCE,
    VIRION_DIAMETER_NM,
    VBOF_REACTION_ID,
    print_configuration,
    # Lipid parameters
    LIPID_FRACTIONS,
    LIPID_PACKING_DENSITY_NM2,
    # Copy-number calculation parameters
    USE_CALCULATED_COPY_NUMBERS,
    VIRION_DIAMETER_NM,
    N_NUCLEOTIDES_PER_PROTOMER,
    GENOME_COPIES_PER_VIRION,
    M_MEMBRANE_COVERAGE_FRACTION,
    M_DIMER_FOOTPRINT_NM2,
    SPIKE_WIDTH_NM,
    SPIKE_SPACING_MODE,
    SPIKE_SPACING_MIN_NM,
    SPIKE_SPACING_AVG_NM,
    SPIKE_SPACING_MAX_NM,
    GLYCOPROTEIN_RATIO,
    F_OLIGOMERIC_STATE,
    G_OLIGOMERIC_STATE,
    SH_OLIGOMERIC_STATE,
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


def main():
    """
    Main pipeline to construct HMPV VBOF.
    
    Steps:
    1. Load and analyze genome
    2. Load and analyze proteins
    3. Calculate stoichiometries
    4. Build VBOF
    5. Save results
    """
    logger.info("=" * 70)
    logger.info("HMPV VBOF CONSTRUCTION PIPELINE")
    logger.info("=" * 70)
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Print configuration
    print_configuration()
    
    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # =====================================================================
        # STEP 1: Load and analyze genome
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 1: GENOME ANALYSIS")
        logger.info("=" * 70)
        
        genome_path = GENOMIC_DIR / DEFAULT_GENOME_FILE
        gff_path = GENOMIC_DIR / DEFAULT_GFF_FILE
        
        # Load genome
        genome = load_genome(str(genome_path))
        genome_summary = get_genome_summary(genome)
        
        # Parse annotations
        annotations = parse_gff_annotations(str(gff_path))
        
        # Calculate genome stoichiometry
        genome_stoichiometry = calculate_genome_stoichiometry(genome, copies_per_virion=1)
        
        logger.info(f"\nGenome Analysis Complete:")
        logger.info(f"  Accession: {genome_summary['accession']}")
        logger.info(f"  Length: {genome_summary['length']} nt")
        logger.info(f"  GC Content: {genome_summary['gc_content']}%")
        logger.info(f"  Genes annotated: {len(annotations)}")
        
        # =====================================================================
        # STEP 2: Load and analyze proteins
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 2: PROTEIN ANALYSIS")
        logger.info("=" * 70)
        
        protein_path = PROTEIN_DIR / DEFAULT_PROTEIN_FILE
        
        # Load proteins
        proteins = load_proteins(str(protein_path))
        protein_summary = get_protein_summary(proteins)
        
        logger.info(f"\nProtein Analysis Complete:")
        logger.info(f"  Total proteins: {protein_summary['total_proteins']}")
        logger.info(f"  Total amino acids: {protein_summary['total_amino_acids']}")
        
        for gene_name, info in protein_summary['proteins'].items():
            logger.info(f"  {gene_name}: {info['length']} aa, MW={info['molecular_weight']} Da")
        
        # =====================================================================
        # STEP 3: Calculate protein stoichiometry with copy numbers
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 3: PROTEIN STOICHIOMETRY")
        logger.info("=" * 70)
        
        # Resolve copy numbers (calculated or default, controlled by USE_CALCULATED_COPY_NUMBERS)
        copy_numbers = get_copy_numbers(
            genome_length=genome.length,
            default_copy_numbers=HMPV_COPY_NUMBERS,
            use_calculated=USE_CALCULATED_COPY_NUMBERS,
            virion_diameter_nm=VIRION_DIAMETER_NM,
            nt_per_protomer=N_NUCLEOTIDES_PER_PROTOMER,
            genome_copies=GENOME_COPIES_PER_VIRION,
            coverage_fraction=M_MEMBRANE_COVERAGE_FRACTION,
            dimer_footprint_nm2=M_DIMER_FOOTPRINT_NM2,
            spike_width_nm=SPIKE_WIDTH_NM,
            spike_spacing_mode=SPIKE_SPACING_MODE,
            spike_spacing_min_nm=SPIKE_SPACING_MIN_NM,
            spike_spacing_avg_nm=SPIKE_SPACING_AVG_NM,
            spike_spacing_max_nm=SPIKE_SPACING_MAX_NM,
            ratio=GLYCOPROTEIN_RATIO,
            oligomeric_states={'F': F_OLIGOMERIC_STATE, 'G': G_OLIGOMERIC_STATE, 'SH': SH_OLIGOMERIC_STATE},
        )

        logger.info(f"\nCopy number source: {'calculated' if USE_CALCULATED_COPY_NUMBERS else 'default'}")
        logger.info("\nUsing copy numbers:")
        for protein, count in copy_numbers.items():
            confidence = COPY_NUMBER_CONFIDENCE.get(protein, 'UNKNOWN')
            logger.info(f"  {protein}: {count} copies/virion (confidence: {confidence})")

        # Calculate protein stoichiometry
        protein_stoichiometry = calculate_protein_stoichiometry(
            proteins,
            copy_numbers=copy_numbers
        )

        # Calculate total amino acids
        total_amino_acids = sum(
            proteins[gene].length * copy_numbers[gene]
            for gene in proteins.keys()
            if gene in copy_numbers
        )

        logger.info(f"\nTotal amino acids in virion: {total_amino_acids}")

        # =====================================================================
        # STEP 4: Build complete VBOF
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 4: VBOF CONSTRUCTION")
        logger.info("=" * 70)

        # Get F and G protein copy numbers for glycan calculation
        f_copy_number = copy_numbers.get('F', 350)
        g_copy_number = copy_numbers.get('G', 250)
        
        vbof = build_vbof(
            genome_stoichiometry=genome_stoichiometry,
            protein_stoichiometry=protein_stoichiometry,
            genome_length=genome.length,
            total_amino_acids=total_amino_acids,
            num_proteins=len(proteins),
            virion_diameter_nm=VIRION_DIAMETER_NM,
            include_lipids=True,
            include_energy=True,
            include_glycans=True,  # NEW: Include glycan requirements
            f_copy_number=f_copy_number,
            g_copy_number=g_copy_number
        )
        
        vbof_summary = get_vbof_summary(vbof)
        
        logger.info(f"\nVBOF Construction Complete:")
        logger.info(f"  Total metabolites: {vbof.metadata['total_metabolites']}")
        logger.info(f"  Consumed: {vbof.metadata['total_consumed']}")
        logger.info(f"  Produced: {vbof.metadata['total_produced']}")
        
        # =====================================================================
        # STEP 5: Save results
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 5: SAVING RESULTS")
        logger.info("=" * 70)
        
        # Save VBOF as JSON
        vbof_dict = export_vbof_to_dict(vbof)
        
        with open(VBOF_JSON_PATH, 'w') as f:
            json.dump(vbof_dict, f, indent=2)
        
        logger.info(f"Saved VBOF to: {VBOF_JSON_PATH}")
        
        # Save summary as text
        with open(VBOF_SUMMARY_PATH, 'w') as f:
            f.write("HMPV VBOF Summary\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("GENOME INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write(f"Accession: {genome_summary['accession']}\n")
            f.write(f"Length: {genome_summary['length']} nt\n")
            f.write(f"GC Content: {genome_summary['gc_content']}%\n\n")
            
            f.write("PROTEIN INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write(f"Copy number source: {'calculated (biophysical)' if USE_CALCULATED_COPY_NUMBERS else 'default (literature)'}\n\n")
            f.write(f"{'Protein':<12} {'Length (aa)':>12} {'Copies/virion':>14} {'Total AAs':>12}\n")
            f.write("-" * 52 + "\n")
            for gene_name, info in protein_summary['proteins'].items():
                copy_num = copy_numbers.get(gene_name, 0)
                total_aas = info['length'] * copy_num if isinstance(copy_num, (int, float)) else 'N/A'
                copy_str = str(copy_num) if copy_num != 'N/A' else 'N/A'
                total_str = f"{total_aas:,}" if isinstance(total_aas, (int, float)) else total_aas
                f.write(f"{gene_name:<12} {info['length']:>12,} {copy_str:>14} {total_str:>12}\n")
            f.write("-" * 52 + "\n")
            f.write(f"{'TOTAL':<12} {'':>12} {'':>14} {total_amino_acids:>12,}\n\n")
            
            # Lipid information
            f.write("LIPID ENVELOPE INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write("Source: Barnes et al. (1987) J Lipid Res 28:130-137 (Sendai virus)\n")
            f.write(f"Virion diameter: {VIRION_DIAMETER_NM} nm\n")
            f.write(f"Lipid packing density: {LIPID_PACKING_DENSITY_NM2} nm²/lipid\n\n")
            lipid_mets = {k: v for k, v in vbof.combined_stoichiometry.items()
                          if k in LIPID_FRACTIONS}
            total_lipids_count = sum(abs(v) for v in lipid_mets.values())
            f.write(f"{'Lipid':<20} {'Fraction':>10} {'Count':>10}\n")
            f.write("-" * 42 + "\n")
            for lipid_id, fraction in LIPID_FRACTIONS.items():
                count = abs(lipid_mets.get(lipid_id, 0))
                f.write(f"{lipid_id:<20} {fraction:>10.4f} {count:>10}\n")
            f.write(f"\nTotal lipid molecules: {total_lipids_count}\n\n")

            # Glycan information
            f.write("GLYCAN INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write("Sources:\n")
            f.write("  F protein: Viswanathan et al. (2011) J Gen Virol 92:1580-1584\n")
            f.write("  G protein: Thammawat et al. (2008) J Virol 82:10022-10034\n\n")
            if 'glycan_data' in vbof.metadata:
                glycan_data = vbof.metadata['glycan_data']
                f.write(f"F protein:\n")
                f.write(f"  N-linked sites: {glycan_data['f_protein']['n_linked_sites']}\n")
                f.write(f"  Total N-glycans: {glycan_data['f_protein']['total_n_glycans']}\n")
                f.write(f"G protein:\n")
                f.write(f"  N-linked sites: {glycan_data['g_protein']['n_linked_sites']}\n")
                f.write(f"  O-linked sites (est.): ~{glycan_data['g_protein']['o_linked_sites_estimated']}\n")
                f.write(f"  Total N-glycans: {glycan_data['g_protein']['total_n_glycans']}\n")
                f.write(f"  Total O-glycans: {glycan_data['g_protein']['total_o_glycans']}\n")
                f.write(f"\nTotal glycan sites: {glycan_data['total_glycan_sites']}\n")
                f.write(f"Glycosylation ATP: {glycan_data['glycosylation_atp']}\n\n")
            
            f.write("VBOF STOICHIOMETRY\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total metabolites: {vbof.metadata['total_metabolites']}\n")
            f.write(f"Consumed: {vbof.metadata['total_consumed']}\n")
            f.write(f"Produced: {vbof.metadata['total_produced']}\n\n")
            
            f.write("CONSUMED METABOLITES (top 20):\n")
            consumed = sorted(
                [(k, v) for k, v in vbof.combined_stoichiometry.items() if v < 0],
                key=lambda x: x[1]
            )
            for met_id, coef in consumed[:20]:
                f.write(f"  {met_id}: {coef}\n")
            
            f.write("\nPRODUCED METABOLITES:\n")
            produced = [(k, v) for k, v in vbof.combined_stoichiometry.items() if v > 0]
            for met_id, coef in produced:
                f.write(f"  {met_id}: {coef}\n")
        
        logger.info(f"Saved summary to: {VBOF_SUMMARY_PATH}")
        
        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 70)
        
        logger.info(f"\nHMPV VBOF successfully constructed!")
       
        
        return vbof
        
    except HMPVModelError as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    vbof = main()

