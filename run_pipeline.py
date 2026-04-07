#!/usr/bin/env python3
"""
HMPV Full Pipeline Orchestrator
================================

Runs the complete HMPV analysis pipeline in sequence:

  1. Build VBOF        — genome + protein stoichiometry → hmpv_vbof.json
  2. Normalize VBOF    — scale to mmol/gDCW → hmpv_vbof_normalized.json
  3. Integrate model   — insert VBOF reaction into host GEM → *_with_HMPV_VBOF.xml
  4. Dual-objective    — gene and reaction knockout analysis (host vs virus)
  5. Selective FVA     — flux variability analysis for shortlisted targets

Single-variant (default)
------------------------
    python run_pipeline.py

  Runs once with VIRION_DIAMETER_NM from config.
  All output goes to OUTPUT_DIR (default: output_HBEC/).

Multi-variant
-------------
    python run_pipeline.py --multi-variant

  Repeats the full pipeline for every entry in VIRION_DIAMETER_VARIANTS (config.py).
  Each variant is isolated under OUTPUT_DIR/variant_<label>/.
  After all variants a cross-variant comparison is written to
  OUTPUT_DIR/multi_variant_comparison/.

Options
-------
    --multi-variant           Run all diameter variants from config
    --variant-diameter FLOAT  Override diameter for a single run (nm)
    --output-dir DIR          Override base output directory
    --virus-max FLOAT         Max virus growth % for selective targets (default 50)
    --host-min FLOAT          Min host growth % for selective targets (default 80)
    --critical FLOAT          Max virus growth % for critical targets (default 5)
    --fraction FLOAT          FVA fraction_of_optimum (default 0.9)
    --no-fva                  Skip step 5
    --no-dual-objective       Skip steps 4 and 5

Author: Syed Mushahid Hussain
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

# --- genome / protein / VBOF construction ---
from src.genome_analyzer import (
    calculate_genome_stoichiometry,
    get_genome_summary,
    load_genome,
    parse_gff_annotations,
)
from src.protein_analyzer import (
    calculate_protein_stoichiometry,
    get_copy_numbers,
    get_protein_summary,
    load_proteins,
)
from src.vbof_builder import build_vbof, export_vbof_to_dict, get_vbof_summary

# --- model integration ---
from src.model_integration import (
    get_integration_summary,
    integrate_vbof,
    load_host_model,
    save_integrated_model,
    validate_integrated_model,
)

# --- dual-objective analysis ---
from src.dual_objective_knockout import (
    DualObjectiveConfig,
    DualObjectiveKnockout,
    load_integrated_model as _load_cobra_model,
)
from src.exceptions import HMPVModelError

# --- config ---
from src.config import (
    # directories and file-name strings
    DUAL_OBJECTIVE_ANALYSIS_DIR,
    MULTI_VARIANT_COMPARISON_DIR,
    # output file names — VBOF / model
    VBOF_JSON_FILE,
    VBOF_NORMALIZED_JSON_FILE,
    VBOF_SUMMARY_FILE,
    INTEGRATION_SUMMARY_FILE,
    INTEGRATED_MODEL_XML_SUFFIX,
    # output file names — dual-objective (gene)
    HOST_GROWTH_RESULTS_FILE,
    VIRUS_GROWTH_RESULTS_FILE,
    MERGED_KNOCKOUT_RESULTS_FILE,
    SELECTIVE_TARGETS_FILE,
    CRITICAL_VIRAL_TARGETS_FILE,
    # output file names — dual-objective (reaction)
    HOST_REACTION_KNOCKOUT_RESULTS_FILE,
    VIRUS_REACTION_KNOCKOUT_RESULTS_FILE,
    MERGED_REACTION_KNOCKOUT_RESULTS_FILE,
    SELECTIVE_REACTION_TARGETS_FILE,
    CRITICAL_REACTION_TARGETS_FILE,
    REACTION_SUBSYSTEM_ESSENTIALITY_FILE,
    # output file names — FVA / comparison
    SELECTIVE_TARGETS_FVA_FILE,
    VARIANT_COMPARISON_SUMMARY_FILE,
    VARIANT_COMPARISON_REPORT_FILE,
    MERGED_SELECTIVE_GENE_FILE,
    MERGED_SELECTIVE_RXN_FILE,
    MERGED_FVA_FILE,
    # input paths
    GENOMIC_DIR,
    PROTEIN_DIR,
    MODEL_DIR,
    OUTPUT_DIR,
    DEFAULT_GENOME_FILE,
    DEFAULT_GFF_FILE,
    DEFAULT_PROTEIN_FILE,
    DEFAULT_HOST_MODEL,
    # reaction IDs
    VBOF_REACTION_ID,
    HOST_BOF_REACTION_ID,
    # VBOF construction parameters
    HMPV_COPY_NUMBERS,
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
    # thresholds
    DEFAULT_THRESHOLDS,
    # variant definitions
    VIRION_DIAMETER_VARIANTS,
)

# --- normalize_vbof functions (root-level script) ---
from normalize_vbof import (
    calculate_virion_mass,
    normalize_vbof as _normalize_coefficients,
)

# --- selective_targets_fva functions (root-level script) ---
from selective_targets_fva import (
    _resolve_biomass_reaction,
    analyze_targets as _run_fva_analysis,
    load_target_tables,
)

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# =============================================================================
# SHARED DATA LOADER  (genome + proteins — invariant across diameter variants)
# =============================================================================

def _load_shared_data() -> tuple:
    """
    Load genome and protein data once; reused for every diameter variant.

    Returns
    -------
    (genome, genome_stoichiometry, genome_summary, proteins, protein_summary)
    """
    genome_path = GENOMIC_DIR / DEFAULT_GENOME_FILE
    gff_path = GENOMIC_DIR / DEFAULT_GFF_FILE
    protein_path = PROTEIN_DIR / DEFAULT_PROTEIN_FILE

    logger.info("Loading genome: %s", genome_path)
    genome = load_genome(str(genome_path))
    parse_gff_annotations(str(gff_path))
    genome_stoich = calculate_genome_stoichiometry(genome, copies_per_virion=1)
    genome_summary = get_genome_summary(genome)

    logger.info("Loading proteins: %s", protein_path)
    proteins = load_proteins(str(protein_path))
    protein_summary = get_protein_summary(proteins)

    logger.info(
        "Genome: %s nt | Proteins loaded: %s",
        genome_summary["length"],
        protein_summary["total_proteins"],
    )
    return genome, genome_stoich, genome_summary, proteins, protein_summary


# =============================================================================
# STEP 1  —  BUILD VBOF
# =============================================================================

def step1_build_vbof(
    genome,
    genome_stoich: dict,
    proteins,
    diameter_nm: float,
    output_dir: Path,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
    genome_copies: int = 1,
) -> dict:
    """
    Build raw VBOF for a given virion geometry.

    Parameters
    ----------
    genome : SeqRecord
        Loaded HMPV genome.
    genome_stoich : dict
        Nucleotide stoichiometry from genome_analyzer.
    proteins : dict
        Loaded HMPV protein sequences.
    diameter_nm : float
        Sphere diameter for spherical; tube cross-section diameter for
        filamentous (nm).  Used for M-protein, glycoprotein, and lipid
        surface-area calculations.
    output_dir : Path
        Variant output directory; hmpv_vbof.json is written here.
    morphology : str
        "spherical" (default) or "filamentous".
        Controls whether 4πr² or 2πrL is used for surface area.
    virion_length_nm : float | None
        Total tip-to-tip length for filamentous particles (nm).
        Required when morphology="filamentous".
    genome_copies : int
        Genome copies per virion (1 for spherical; 2 for filamentous).

    Returns
    -------
    dict  Raw VBOF data (same structure as export_vbof_to_dict output).
    """
    logger.info("\n" + "=" * 70)
    logger.info(
        "STEP 1: BUILD VBOF  (%s  diameter=%.1f nm%s)",
        morphology,
        diameter_nm,
        f"  length={virion_length_nm:.1f} nm" if virion_length_nm else "",
    )
    logger.info("=" * 70)

    copy_numbers = get_copy_numbers(
        genome_length=genome.length,
        default_copy_numbers=HMPV_COPY_NUMBERS,
        use_calculated=USE_CALCULATED_COPY_NUMBERS,
        virion_diameter_nm=diameter_nm,
        nt_per_protomer=N_NUCLEOTIDES_PER_PROTOMER,
        genome_copies=genome_copies,
        coverage_fraction=M_MEMBRANE_COVERAGE_FRACTION,
        dimer_footprint_nm2=M_DIMER_FOOTPRINT_NM2,
        spike_width_nm=SPIKE_WIDTH_NM,
        spike_spacing_mode=SPIKE_SPACING_MODE,
        spike_spacing_min_nm=SPIKE_SPACING_MIN_NM,
        spike_spacing_avg_nm=SPIKE_SPACING_AVG_NM,
        spike_spacing_max_nm=SPIKE_SPACING_MAX_NM,
        ratio=GLYCOPROTEIN_RATIO,
        oligomeric_states={
            "F": F_OLIGOMERIC_STATE,
            "G": G_OLIGOMERIC_STATE,
            "SH": SH_OLIGOMERIC_STATE,
        },
        morphology=morphology,
        virion_length_nm=virion_length_nm,
    )

    protein_stoich = calculate_protein_stoichiometry(proteins, copy_numbers=copy_numbers)
    total_aa = sum(
        proteins[g].length * copy_numbers[g]
        for g in proteins
        if g in copy_numbers
    )

    vbof = build_vbof(
        genome_stoichiometry=genome_stoich,
        protein_stoichiometry=protein_stoich,
        genome_length=genome.length,
        total_amino_acids=total_aa,
        num_proteins=len(proteins),
        virion_diameter_nm=diameter_nm,
        include_lipids=True,
        include_energy=True,
        include_glycans=True,
        f_copy_number=copy_numbers.get("F", 695),
        g_copy_number=copy_numbers.get("G", 139),
        morphology=morphology,
        virion_length_nm=virion_length_nm,
    )

    vbof_dict = export_vbof_to_dict(vbof)
    vbof_path = output_dir / VBOF_JSON_FILE
    with open(vbof_path, "w") as fh:
        json.dump(vbof_dict, fh, indent=2)

    logger.info(
        "Saved VBOF → %s  (%s metabolites)",
        vbof_path,
        vbof.metadata.get("total_metabolites", "?"),
    )

    # Write human-readable summary
    summary = get_vbof_summary(vbof)
    summary_path = output_dir / VBOF_SUMMARY_FILE
    with open(summary_path, "w") as fh:
        fh.write(f"HMPV VBOF Summary\n{'=' * 70}\n\n")
        fh.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"Morphology: {morphology}  |  Diameter: {diameter_nm} nm")
        if virion_length_nm:
            fh.write(f"  |  Length: {virion_length_nm} nm")
        fh.write(f"\nGenome copies per virion: {genome_copies}\n\n")
        meta = summary["metadata"]
        fh.write(f"Total metabolites : {meta.get('total_metabolites', '?')}\n")
        fh.write(f"Consumed          : {meta.get('total_consumed', '?')}\n")
        fh.write(f"Produced          : {meta.get('total_produced', '?')}\n\n")
        for cat, info in summary["categories"].items():
            fh.write(f"  {cat:<14s}: {info['count']} metabolites  "
                     f"(total coef {info['total']:.4e})\n")
        fh.write("\nCONSUMED METABOLITES (top 20):\n")
        for met, coef in sorted(summary["consumed"].items(), key=lambda x: x[1])[:20]:
            fh.write(f"  {met}: {coef}\n")
        fh.write("\nPRODUCED METABOLITES:\n")
        for met, coef in summary["produced"].items():
            fh.write(f"  {met}: {coef}\n")
    logger.info("Saved VBOF summary → %s", summary_path)

    return vbof_dict


# =============================================================================
# STEP 2  —  NORMALIZE VBOF
# =============================================================================

def step2_normalize_vbof(vbof_dict: dict, output_dir: Path) -> dict:
    """
    Normalize raw VBOF coefficients to mmol / gDCW of virion biomass.

    Parameters
    ----------
    vbof_dict : dict
        Raw VBOF dict as returned by step1_build_vbof.
    output_dir : Path
        Variant output directory; hmpv_vbof_normalized.json is written here.

    Returns
    -------
    dict  Normalized VBOF dict (key 'combined_stoichiometry' holds mmol/gDCW values).
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: NORMALIZE VBOF")
    logger.info("=" * 70)

    raw_stoich = vbof_dict["combined_stoichiometry"]
    total_consumed = sum(abs(c) for c in raw_stoich.values() if c < 0)
    virion_mass = calculate_virion_mass(raw_stoich)
    normalized = _normalize_coefficients(raw_stoich, target_mass_g=1.0)

    normalized_dict = {
        "combined_stoichiometry": normalized,
        "metadata": {
            **vbof_dict.get("metadata", {}),
            "normalization_method": "virion_mass_mmol_per_gDCW",
            "normalization_factor": total_consumed,
            "virion_mass_g": virion_mass,
            "normalized_at": datetime.now().isoformat(),
        },
        "raw_genome_stoichiometry": vbof_dict.get("genome_stoichiometry", {}),
        "raw_protein_stoichiometry": vbof_dict.get("protein_stoichiometry", {}),
        "raw_energy_stoichiometry": vbof_dict.get("energy_stoichiometry", {}),
        "raw_lipid_stoichiometry": vbof_dict.get("lipid_stoichiometry", {}),
        "raw_glycan_stoichiometry": vbof_dict.get("glycan_stoichiometry", {}),
    }

    norm_path = output_dir / VBOF_NORMALIZED_JSON_FILE
    with open(norm_path, "w") as fh:
        json.dump(normalized_dict, fh, indent=2)

    logger.info(
        "Saved normalized VBOF → %s  (virion mass = %.4e g)", norm_path, virion_mass
    )
    return normalized_dict


# =============================================================================
# STEP 3  —  INTEGRATE MODEL
# =============================================================================

def step3_integrate_model(normalized_dict: dict, output_dir: Path) -> Path:
    """
    Integrate the normalized VBOF reaction into the host GEM.

    Parameters
    ----------
    normalized_dict : dict
        Normalized VBOF dict as returned by step 2 normalize_vbof.
    output_dir : Path
        Variant output directory; {model_id}_with_HMPV_VBOF.xml is written here.

    Returns
    -------
    Path  Absolute path to the saved integrated model XML.
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: INTEGRATE MODEL")
    logger.info("=" * 70)

    stoich = normalized_dict["combined_stoichiometry"]
    host_model = load_host_model(MODEL_DIR / DEFAULT_HOST_MODEL)
    logger.info(
        "Host model: %s  (%s reactions, %s genes)",
        host_model.id, len(host_model.reactions), len(host_model.genes),
    )

    integrated_model, _vbof_rxn, unmapped = integrate_vbof(
        model=host_model,
        vbof_stoichiometry=stoich,
        reaction_id=VBOF_REACTION_ID,
        skip_unmapped=True,
    )

    if unmapped:
        logger.warning(
            "%s unmapped metabolites skipped (first 5: %s)",
            len(unmapped), list(unmapped)[:5],
        )

    validation = validate_integrated_model(integrated_model, VBOF_REACTION_ID)

    model_xml = output_dir / f"{integrated_model.id}{INTEGRATED_MODEL_XML_SUFFIX}"
    save_integrated_model(integrated_model, model_xml, format="sbml")
    logger.info("Saved integrated model → %s", model_xml)

    # Write integration summary
    summary_text = get_integration_summary(
        model=integrated_model,
        vbof_reaction_id=VBOF_REACTION_ID,
        unmapped=list(unmapped),
        validation=validation,
    )
    summary_path = output_dir / INTEGRATION_SUMMARY_FILE
    with open(summary_path, "w") as fh:
        fh.write(summary_text)
    logger.info("Saved integration summary → %s", summary_path)

    return model_xml


# =============================================================================
# STEP 4  —  DUAL-OBJECTIVE KNOCKOUT ANALYSIS
# =============================================================================

def step4_dual_objective(
    model_xml: Path,
    output_dir: Path,
    virus_max: float,
    host_min: float,
    critical_pct: float,
) -> dict:
    """
    Run dual-objective gene and reaction knockout analysis.

    Parameters
    ----------
    model_xml : Path
        Integrated model XML (output of step3).
    output_dir : Path
        Variant root directory; results are saved under
        output_dir/dual_objective_analysis/.
    virus_max : float
        Max virus growth % threshold for selective target classification.
    host_min : float
        Min host growth % threshold for selective target classification.
    critical_pct : float
        Max virus growth % threshold for critical target classification.

    Returns
    -------
    dict  Metadata including result file paths and summary counts.
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 4: DUAL-OBJECTIVE KNOCKOUT ANALYSIS")
    logger.info(
        "  thresholds: virus < %.0f%%  |  host > %.0f%%  |  critical < %.0f%%",
        virus_max, host_min, critical_pct,
    )
    logger.info("=" * 70)

    duo_dir = output_dir / DUAL_OBJECTIVE_ANALYSIS_DIR
    duo_dir.mkdir(parents=True, exist_ok=True)

    model = _load_cobra_model(model_xml, vbof_id=VBOF_REACTION_ID)

    thresholds = DEFAULT_THRESHOLDS.copy()
    thresholds["selective_target"] = {
        "virus_max": virus_max / 100.0,
        "host_min": host_min / 100.0,
    }
    thresholds["lethal_virus"] = {"virus_max": critical_pct / 100.0}

    cfg = DualObjectiveConfig(
        vbof_id=VBOF_REACTION_ID,
        bof_id=HOST_BOF_REACTION_ID,
        thresholds=thresholds,
    )
    analyzer = DualObjectiveKnockout(model, cfg)

    # --- gene-level knockouts ---
    host_results = analyzer.run_host_knockout_analysis()
    virus_results = analyzer.run_virus_knockout_analysis()
    merged = analyzer.merge_results(host_results, virus_results)
    selective = analyzer.filter_selective_targets(merged, virus_max=virus_max, host_min=host_min)
    critical = analyzer.filter_critical_viral_targets(merged, virus_max=critical_pct)

    # --- reaction-level knockouts ---
    host_rxn = analyzer.run_host_reaction_knockout_analysis()
    virus_rxn = analyzer.run_virus_reaction_knockout_analysis()
    merged_rxn = analyzer.merge_reaction_results(host_rxn, virus_rxn)
    selective_rxn = analyzer.filter_selective_reaction_targets(
        merged_rxn, virus_max=virus_max, host_min=host_min
    )
    critical_rxn = analyzer.filter_critical_reaction_targets(merged_rxn, virus_max=critical_pct)
    subsystem = analyzer.analyze_subsystem_essentiality(merged_rxn)

    # --- save all CSVs ---
    def _save(df: pd.DataFrame, fname: str) -> Path:
        p = duo_dir / fname
        df.to_csv(p, index=False)
        logger.info("  %-50s  (%s rows)", fname, len(df))
        return p

    _save(host_results, HOST_GROWTH_RESULTS_FILE)
    _save(virus_results, VIRUS_GROWTH_RESULTS_FILE)
    _save(merged, MERGED_KNOCKOUT_RESULTS_FILE)
    _save(selective, SELECTIVE_TARGETS_FILE)
    _save(critical, CRITICAL_VIRAL_TARGETS_FILE)
    _save(host_rxn, HOST_REACTION_KNOCKOUT_RESULTS_FILE)
    _save(virus_rxn, VIRUS_REACTION_KNOCKOUT_RESULTS_FILE)
    _save(merged_rxn, MERGED_REACTION_KNOCKOUT_RESULTS_FILE)
    _save(selective_rxn, SELECTIVE_REACTION_TARGETS_FILE)
    _save(critical_rxn, CRITICAL_REACTION_TARGETS_FILE)
    _save(subsystem, REACTION_SUBSYSTEM_ESSENTIALITY_FILE)

    logger.info(
        "Dual-objective complete: %s selective gene targets | %s selective reaction targets",
        len(selective), len(selective_rxn),
    )

    return {
        "duo_dir": duo_dir,
        "selective_targets_path": duo_dir / SELECTIVE_TARGETS_FILE,
        "selective_rxn_targets_path": duo_dir / SELECTIVE_REACTION_TARGETS_FILE,
        "n_selective_gene": len(selective),
        "n_selective_rxn": len(selective_rxn),
        "host_baseline": getattr(analyzer, "host_baseline", float("nan")),
        "virus_baseline": getattr(analyzer, "virus_baseline", float("nan")),
        "n_reactions": len(model.reactions),
        "n_genes": len(model.genes),
    }


# =============================================================================
# STEP 5  —  SELECTIVE TARGETS FVA
# =============================================================================

def step5_selective_fva(
    model_xml: Path,
    selective_rxn_path: Path,
    selective_gene_path: Path,
    output_dir: Path,
    fraction: float = 0.9,
    run_host: bool = True,
) -> Optional[Path]:
    """
    Run flux variability analysis on the shortlisted selective targets.

    Parameters
    ----------
    model_xml : Path
        Integrated model XML (fresh load — no accumulated knockouts).
    selective_rxn_path : Path
        selective_reaction_targets.csv from step4.
    selective_gene_path : Path
        selective_antiviral_targets.csv from step4.
    output_dir : Path
        Variant root directory; FVA CSV is saved under
        output_dir/dual_objective_analysis/.
    fraction : float
        fraction_of_optimum for COBRApy FVA (default 0.9).
    run_host : bool
        If True, also run FVA with host BOF objective.

    Returns
    -------
    Path | None  Path to the saved FVA CSV, or None if no targets were loaded.
    """
    logger.info("\n" + "=" * 70)
    logger.info(
        "STEP 5: SELECTIVE TARGETS FVA  (fraction=%.2f | host FVA=%s)",
        fraction, "on" if run_host else "off",
    )
    logger.info("=" * 70)

    model = _load_cobra_model(model_xml, vbof_id=VBOF_REACTION_ID)

    if VBOF_REACTION_ID not in model.reactions:
        raise KeyError(f"VBOF reaction {VBOF_REACTION_ID!r} not in model.")
    _resolve_biomass_reaction(model, HOST_BOF_REACTION_ID)

    targets = load_target_tables(selective_rxn_path, selective_gene_path)
    if targets.empty:
        logger.warning("No FVA targets found — skipping FVA step.")
        return None

    logger.info(
        "Running FVA for %s targets (%s reaction, %s gene) …",
        len(targets),
        (targets["knockout_type"] == "reaction").sum(),
        (targets["knockout_type"] == "gene").sum(),
    )

    result = _run_fva_analysis(
        model,
        targets,
        vbof_id=VBOF_REACTION_ID,
        host_bof_id=HOST_BOF_REACTION_ID,
        fraction=fraction,
        run_host=run_host,
    )

    col_order = [
        "target_id", "knockout_type",
        "virus_fva_min", "virus_fva_max",
        "host_fva_min", "host_fva_max",
        "status",
    ]
    result = result[[c for c in col_order if c in result.columns]]

    duo_dir = output_dir / DUAL_OBJECTIVE_ANALYSIS_DIR
    duo_dir.mkdir(parents=True, exist_ok=True)
    fva_path = duo_dir / SELECTIVE_TARGETS_FVA_FILE
    result.to_csv(fva_path, index=False)

    ok_count = (result["status"] == "ok").sum()
    logger.info(
        "FVA complete: %s/%s successful → %s", ok_count, len(result), fva_path
    )
    return fva_path


# =============================================================================
# SINGLE-VARIANT RUNNER
# =============================================================================

def _run_one_variant(
    args: argparse.Namespace,
    genome,
    genome_stoich: dict,
    proteins,
    diameter_nm: float,
    label: str,
    output_dir: Path,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
    genome_copies: int = 1,
) -> dict:
    """
    Run the complete 5-step pipeline for one virion geometry variant.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (thresholds, flags).
    genome, genome_stoich, proteins : shared data from _load_shared_data().
    diameter_nm : float
        Sphere diameter (spherical) or tube cross-section diameter (filamentous) in nm.
    label : str
        Short identifier used in logging and output paths.
    output_dir : Path
        Variant-specific output directory.
    morphology : str
        "spherical" or "filamentous" — selects the surface-area formula.
    virion_length_nm : float | None
        Total tip-to-tip length for filamentous particles (nm).
    genome_copies : int
        Genome copies per virion (1 for spherical, 2 for filamentous).

    Returns
    -------
    dict  Per-variant summary metrics for the comparison report.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("\n" + "#" * 70)
    logger.info(
        "PIPELINE  variant=%s  morphology=%s  diameter=%.1f nm%s",
        label, morphology, diameter_nm,
        f"  length={virion_length_nm:.1f} nm" if virion_length_nm else "",
    )
    logger.info("#" * 70)

    # Step 1 — build  (pass morphology + length so correct surface-area formula is used)
    vbof_dict = step1_build_vbof(
        genome, genome_stoich, proteins,
        diameter_nm=diameter_nm,
        output_dir=output_dir,
        morphology=morphology,
        virion_length_nm=virion_length_nm,
        genome_copies=genome_copies,
    )

    # Step 2 — normalize
    norm_dict = step2_normalize_vbof(vbof_dict, output_dir)

    # Step 3 — integrate
    model_xml = step3_integrate_model(norm_dict, output_dir)

    summary: dict = {
        "variant": label,
        "morphology": morphology,
        "diameter_nm": diameter_nm,
        "length_nm": virion_length_nm if virion_length_nm else "N/A",
        "genome_copies": genome_copies,
        "vbof_metabolites": vbof_dict.get("metadata", {}).get("total_metabolites", "?"),
        "n_selective_gene": 0,
        "n_selective_rxn": 0,
        "n_fva_ok": 0,
        "host_baseline": float("nan"),
        "virus_baseline": float("nan"),
        "n_reactions": "?",
        "n_genes": "?",
        "output_dir": str(output_dir),
        "status": "ok",
    }

    if args.no_dual_objective:
        logger.info("Skipping steps 4 & 5 (--no-dual-objective).")
        return summary

    # Step 4 — dual-objective
    duo_info = step4_dual_objective(
        model_xml,
        output_dir,
        virus_max=args.virus_max,
        host_min=args.host_min,
        critical_pct=args.critical,
    )
    summary.update({
        "n_selective_gene": duo_info["n_selective_gene"],
        "n_selective_rxn": duo_info["n_selective_rxn"],
        "host_baseline": duo_info["host_baseline"],
        "virus_baseline": duo_info["virus_baseline"],
        "n_reactions": duo_info["n_reactions"],
        "n_genes": duo_info["n_genes"],
    })

    if args.no_fva:
        logger.info("Skipping step 5 (--no-fva).")
        return summary

    # Step 5 — FVA
    fva_path = step5_selective_fva(
        model_xml,
        selective_rxn_path=duo_info["selective_rxn_targets_path"],
        selective_gene_path=duo_info["selective_targets_path"],
        output_dir=output_dir,
        fraction=args.fraction,
        run_host=True,
    )
    if fva_path and fva_path.exists():
        fva_df = pd.read_csv(fva_path)
        summary["n_fva_ok"] = int((fva_df["status"] == "ok").sum())

    return summary


# =============================================================================
# MULTI-VARIANT ORCHESTRATOR
# =============================================================================

def _collect_variant_csv(
    variant_dir: Path,
    subdirectory: str,
    filename: str,
    label: str,
    diameter_nm: float,
    length_nm: Optional[float] = None,
) -> Optional[pd.DataFrame]:
    """
    Load a variant result CSV and prepend variant/diameter/length columns.
    length_nm is NaN for spherical particles (not applicable).
    Returns None if the file does not exist.
    """
    path = variant_dir / subdirectory / filename
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df.insert(0, "variant", label)
    df.insert(1, "diameter_nm", diameter_nm)
    len_val = float("nan") if length_nm is None else float(length_nm)
    df.insert(2, "length_nm", len_val)
    return df


def _write_comparison_report(
    summaries: List[dict],
    merged_gene: Optional[pd.DataFrame],
    merged_rxn: Optional[pd.DataFrame],
    args: argparse.Namespace,
    comp_dir: Path,
) -> None:
    """Write the cross-variant text comparison report."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_variants = len(summaries)

    lines: List[str] = [
        "=" * 80,
        "HMPV MULTI-VARIANT COMPARISON REPORT",
        "=" * 80,
        f"Generated  : {now}",
        f"Variants   : {n_variants}",
        f"Thresholds : virus < {args.virus_max}%  |  "
        f"host > {args.host_min}%  |  critical < {args.critical}%",
        "",
        "=" * 80,
        "VARIANT OVERVIEW",
        "=" * 80,
        f"{'Variant':<28} {'Morphology':<13} {'Diam(nm)':>9} {'Len(nm)':>8} "
        f"{'GeneTgts':>9} {'RxnTgts':>8} {'FVA OK':>7} {'HostFlux':>12} {'VirusFlux':>12}",
        "-" * 110,
    ]

    for s in summaries:
        host_flux  = s.get("host_baseline", float("nan"))
        virus_flux = s.get("virus_baseline", float("nan"))
        length_val = s.get("length_nm", "N/A")
        length_str = f"{length_val:.0f}" if isinstance(length_val, float) else str(length_val)
        lines.append(
            f"{s.get('variant', '?'):<28} "
            f"{s.get('morphology', '?'):<13} "
            f"{s.get('diameter_nm', 0):>9.1f} "
            f"{length_str:>8} "
            f"{str(s.get('n_selective_gene', '-')):>9} "
            f"{str(s.get('n_selective_rxn', '-')):>8} "
            f"{str(s.get('n_fva_ok', '-')):>7} "
            f"{host_flux:>12.6f} "
            f"{virus_flux:>12.6f}"
        )

    # --- conserved gene targets ---
    if merged_gene is not None and "gene_id" in merged_gene.columns:
        gene_counts = merged_gene.groupby("gene_id")["variant"].nunique()
        conserved = gene_counts[gene_counts == n_variants].index.tolist()
        lines += [
            "",
            "=" * 80,
            f"CONSERVED SELECTIVE GENE TARGETS  (present in all {n_variants} variants)",
            "=" * 80,
        ]
        if conserved:
            for g in conserved[:40]:
                lines.append(f"  {g}")
            if len(conserved) > 40:
                lines.append(f"  ... and {len(conserved) - 40} more")
        else:
            lines.append("  None shared across all variants.")

        # unique per variant
        lines += ["", "Targets unique to each variant:", ""]
        for s in summaries:
            lbl = s.get("variant", "?")
            diam = s.get("diameter_nm", "?")
            variant_genes = set(
                merged_gene.loc[merged_gene["variant"] == lbl, "gene_id"]
            )
            unique = [
                g for g in variant_genes
                if gene_counts.get(g, 0) == 1
            ]
            lines.append(f"  {lbl} ({diam} nm) — {len(unique)} unique targets")
            for g in unique[:10]:
                lines.append(f"    {g}")

    # --- conserved reaction targets ---
    if merged_rxn is not None and "reaction_id" in merged_rxn.columns:
        rxn_counts = merged_rxn.groupby("reaction_id")["variant"].nunique()
        conserved_rxn = rxn_counts[rxn_counts == n_variants].index.tolist()
        lines += [
            "",
            "=" * 80,
            f"CONSERVED SELECTIVE REACTION TARGETS  (present in all {n_variants} variants)",
            "=" * 80,
        ]
        if conserved_rxn:
            for r in conserved_rxn[:40]:
                lines.append(f"  {r}")
            if len(conserved_rxn) > 40:
                lines.append(f"  ... and {len(conserved_rxn) - 40} more")
        else:
            lines.append("  None shared across all variants.")

    # --- output files summary ---
    lines += [
        "",
        "=" * 80,
        "OUTPUT FILES IN THIS DIRECTORY",
        "=" * 80,
        f"  {VARIANT_COMPARISON_SUMMARY_FILE:<45} per-variant metric table",
        f"  {VARIANT_COMPARISON_REPORT_FILE:<45} this report",
        f"  {MERGED_SELECTIVE_GENE_FILE:<45} all selective gene targets (variant column added)",
        f"  {MERGED_SELECTIVE_RXN_FILE:<45} all selective reaction targets (variant column added)",
        f"  {MERGED_FVA_FILE:<45} all FVA results (variant column added)",
        "",
        "Per-variant results are in:",
    ]
    for s in summaries:
        lines.append(f"  {s.get('output_dir', '?')}")

    lines += ["", "=" * 80]

    report_path = comp_dir / VARIANT_COMPARISON_REPORT_FILE
    report_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("Saved comparison report → %s", report_path)


def run_multi_variant(
    args: argparse.Namespace,
    genome,
    genome_stoich: dict,
    proteins,
    base_output_dir: Path,
) -> None:
    """
    Iterate over all VIRION_DIAMETER_VARIANTS, run the complete pipeline for
    each one, then collect results and write a cross-variant comparison.
    """
    variants: Dict[str, Dict] = VIRION_DIAMETER_VARIANTS
    if not variants:
        logger.error(
            "VIRION_DIAMETER_VARIANTS is empty in config.py — add at least one variant."
        )
        sys.exit(1)

    logger.info("\n" + "#" * 70)
    logger.info("MULTI-VARIANT MODE  —  %s variants", len(variants))
    for lbl, vp in variants.items():
        d = vp["diameter_nm"]
        ln = vp.get("length_nm")
        extra = f"  L={float(ln):.0f} nm" if ln is not None else ""
        logger.info("  %-30s Ø %.1f nm%s", lbl, d, extra)
    logger.info("#" * 70)

    all_summaries: List[dict] = []
    all_gene_frames: List[pd.DataFrame] = []
    all_rxn_frames: List[pd.DataFrame] = []
    all_fva_frames: List[pd.DataFrame] = []

    for label, vparams in variants.items():
        diameter_nm    = vparams["diameter_nm"]
        morphology     = vparams.get("morphology", "spherical")
        virion_length_nm = vparams.get("length_nm")
        genome_copies  = vparams.get("genome_copies", 1)

        variant_dir = base_output_dir / f"variant_{label}"
        try:
            summary = _run_one_variant(
                args, genome, genome_stoich, proteins,
                diameter_nm=diameter_nm,
                label=label,
                output_dir=variant_dir,
                morphology=morphology,
                virion_length_nm=virion_length_nm,
                genome_copies=genome_copies,
            )
            all_summaries.append(summary)

            # collect gene targets
            df = _collect_variant_csv(
                variant_dir, DUAL_OBJECTIVE_ANALYSIS_DIR,
                SELECTIVE_TARGETS_FILE, label, diameter_nm,
                length_nm=virion_length_nm,
            )
            if df is not None:
                all_gene_frames.append(df)

            # collect reaction targets
            df = _collect_variant_csv(
                variant_dir, DUAL_OBJECTIVE_ANALYSIS_DIR,
                SELECTIVE_REACTION_TARGETS_FILE, label, diameter_nm,
                length_nm=virion_length_nm,
            )
            if df is not None:
                all_rxn_frames.append(df)

            # collect FVA results
            df = _collect_variant_csv(
                variant_dir, DUAL_OBJECTIVE_ANALYSIS_DIR,
                SELECTIVE_TARGETS_FVA_FILE, label, diameter_nm,
                length_nm=virion_length_nm,
            )
            if df is not None:
                all_fva_frames.append(df)

        except Exception as exc:
            logger.exception("Variant '%s' FAILED: %s", label, exc)
            all_summaries.append({
                "variant": label,
                "morphology": morphology,
                "diameter_nm": diameter_nm,
                "length_nm": virion_length_nm if virion_length_nm else "N/A",
                "status": f"ERROR: {exc}",
                "output_dir": str(variant_dir),
            })

    # --- write comparison outputs ---
    comp_dir = base_output_dir / MULTI_VARIANT_COMPARISON_DIR
    comp_dir.mkdir(parents=True, exist_ok=True)

    logger.info("\n" + "=" * 70)
    logger.info("WRITING MULTI-VARIANT COMPARISON  →  %s", comp_dir)
    logger.info("=" * 70)

    # Summary CSV (one row per variant)
    pd.DataFrame(all_summaries).to_csv(
        comp_dir / VARIANT_COMPARISON_SUMMARY_FILE, index=False
    )
    logger.info("  %s", VARIANT_COMPARISON_SUMMARY_FILE)

    merged_gene = None
    merged_rxn = None

    # Merged gene targets
    if all_gene_frames:
        merged_gene = pd.concat(all_gene_frames, ignore_index=True)
        merged_gene.to_csv(comp_dir / MERGED_SELECTIVE_GENE_FILE, index=False)
        logger.info("  %s  (%s total rows)", MERGED_SELECTIVE_GENE_FILE, len(merged_gene))

    # Merged reaction targets
    if all_rxn_frames:
        merged_rxn = pd.concat(all_rxn_frames, ignore_index=True)
        merged_rxn.to_csv(comp_dir / MERGED_SELECTIVE_RXN_FILE, index=False)
        logger.info("  %s  (%s total rows)", MERGED_SELECTIVE_RXN_FILE, len(merged_rxn))

    # Merged FVA results
    if all_fva_frames:
        merged_fva = pd.concat(all_fva_frames, ignore_index=True)
        merged_fva.to_csv(comp_dir / MERGED_FVA_FILE, index=False)
        logger.info("  %s  (%s total rows)", MERGED_FVA_FILE, len(merged_fva))

    # Text comparison report
    _write_comparison_report(all_summaries, merged_gene, merged_rxn, args, comp_dir)

    logger.info("\nMulti-variant analysis complete.")
    logger.info("Results in: %s", comp_dir)


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- mode ---
    p.add_argument(
        "--multi-variant", action="store_true",
        help=(
            "Run the full pipeline for every diameter variant in "
            "VIRION_DIAMETER_VARIANTS (config.py) and generate a "
            "cross-variant comparison report."
        ),
    )
    p.add_argument(
        "--variant-diameter", type=float, default=None, metavar="NM",
        help=(
            "Virion diameter in nm for a single run "
            f"(default: VIRION_DIAMETER_NM = {VIRION_DIAMETER_NM} nm from config)."
        ),
    )

    # --- paths ---
    p.add_argument(
        "--output-dir", type=Path, default=OUTPUT_DIR, metavar="DIR",
        help=f"Base output directory (default: {OUTPUT_DIR}).",
    )

    # --- selectivity thresholds ---
    p.add_argument(
        "--virus-max", type=float, default=50.0, metavar="PCT",
        help="Max virus growth %% for selective targets (default: 50).",
    )
    p.add_argument(
        "--host-min", type=float, default=80.0, metavar="PCT",
        help="Min host growth %% for selective targets (default: 80).",
    )
    p.add_argument(
        "--critical", type=float, default=5.0, metavar="PCT",
        help="Max virus growth %% for critical targets (default: 5).",
    )

    # --- FVA ---
    p.add_argument(
        "--fraction", type=float, default=0.9, metavar="F",
        help="FVA fraction_of_optimum (default: 0.9).",
    )

    # --- skip flags ---
    p.add_argument(
        "--no-fva", action="store_true",
        help="Skip step 5 (selective targets FVA).",
    )
    p.add_argument(
        "--no-dual-objective", action="store_true",
        help="Skip steps 4 and 5 (dual-objective analysis and FVA).",
    )

    return p.parse_args()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    args = parse_args()

    logger.info("=" * 70)
    logger.info("HMPV FULL PIPELINE ORCHESTRATOR")
    logger.info("Started : %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("=" * 70)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Genome + proteins are diameter-independent — load once and share
    logger.info("\n" + "=" * 70)
    logger.info("LOADING SHARED GENOME + PROTEIN DATA")
    logger.info("=" * 70)
    genome, genome_stoich, _genome_summary, proteins, _protein_summary = (
        _load_shared_data()
    )

    if args.multi_variant:
        run_multi_variant(args, genome, genome_stoich, proteins, args.output_dir)
    else:
        diameter_nm = (
            args.variant_diameter
            if args.variant_diameter is not None
            else VIRION_DIAMETER_NM
        )
        label = f"spherical_{diameter_nm:.0f}nm"
        _run_one_variant(
            args, genome, genome_stoich, proteins,
            diameter_nm=diameter_nm,
            label=label,
            output_dir=args.output_dir,
            morphology="spherical",
            virion_length_nm=None,
            genome_copies=GENOME_COPIES_PER_VIRION,
        )

    logger.info("\n" + "=" * 70)
    logger.info("PIPELINE COMPLETE")
    logger.info("Finished: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
