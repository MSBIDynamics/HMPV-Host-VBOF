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
    get_protein_summary
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
    HMPV_COPY_NUMBERS,
    COPY_NUMBER_CONFIDENCE,
    VIRION_DIAMETER_NM,
    VBOF_REACTION_ID,
    print_configuration
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
        
        logger.info("\nUsing copy numbers:")
        for protein, count in HMPV_COPY_NUMBERS.items():
            confidence = COPY_NUMBER_CONFIDENCE.get(protein, 'UNKNOWN')
            logger.info(f"  {protein}: {count} copies/virion (confidence: {confidence})")
        
        # Calculate protein stoichiometry
        protein_stoichiometry = calculate_protein_stoichiometry(
            proteins,
            copy_numbers=HMPV_COPY_NUMBERS
        )
        
        # Calculate total amino acids
        total_amino_acids = sum(
            proteins[gene].length * HMPV_COPY_NUMBERS[gene]
            for gene in proteins.keys()
            if gene in HMPV_COPY_NUMBERS
        )
        
        logger.info(f"\nTotal amino acids in virion: {total_amino_acids}")
        
        # =====================================================================
        # STEP 4: Build complete VBOF
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("STEP 4: VBOF CONSTRUCTION")
        logger.info("=" * 70)
        
        # Get F and G protein copy numbers for glycan calculation
        f_copy_number = HMPV_COPY_NUMBERS.get('F', 350)
        g_copy_number = HMPV_COPY_NUMBERS.get('G', 250)
        
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
        vbof_json_path = output_dir / "hmpv_vbof.json"
        
        with open(vbof_json_path, 'w') as f:
            json.dump(vbof_dict, f, indent=2)
        
        logger.info(f"Saved VBOF to: {vbof_json_path}")
        
        # Save summary as text
        summary_path = output_dir / "hmpv_vbof_summary.txt"
        with open(summary_path, 'w') as f:
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
            for gene_name, info in protein_summary['proteins'].items():
                copy_num = HMPV_COPY_NUMBERS.get(gene_name, 'N/A')
                f.write(f"{gene_name}: {info['length']} aa, {copy_num} copies/virion\n")
            f.write(f"\nTotal amino acids in virion: {total_amino_acids}\n\n")
            
            # NEW: Add glycan information
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
        
        logger.info(f"Saved summary to: {summary_path}")
        
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

