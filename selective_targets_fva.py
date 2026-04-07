#!/usr/bin/env python3
"""
Flux Variability Analysis (FVA) for Selective Dual-Objective Targets
=====================================================================

Reads reaction-level and gene-level selective target tables produced by
dual_objective_analysis.py, applies each knockout on the integrated host+VBOF
model, and runs COBRApy FVA (fraction_of_optimum=0.9) with the viral biomass
reaction (HMPV_VBOF) and optionally the host biomass reaction as objectives.

Each knockout is applied inside ``with model:`` so the base model is restored
automatically after every target.

Output Files:
-------------
- selective_targets_fva.csv: FVA flux ranges for each selective target

Usage:
------
    python selective_targets_fva.py
    python selective_targets_fva.py --model output_HBEC/model_with_HMPV_VBOF.xml
    python selective_targets_fva.py --no-host-fva

Author: Syed Mushahid Hussain
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from cobra import Model
from cobra.flux_analysis import flux_variability_analysis

# Project imports (script lives at repo root)
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.config import (
    HOST_BOF_REACTION_ID,
    INTEGRATED_MODEL_XML_SUFFIX,
    OUTPUT_DIR,
    SELECTIVE_REACTION_TARGETS_PATH,
    SELECTIVE_TARGETS_FVA_PATH,
    SELECTIVE_TARGETS_PATH,
    VBOF_REACTION_ID,
)
from src.dual_objective_knockout import load_integrated_model

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def resolve_integrated_model_path(explicit: Optional[Path]) -> Path:
    """
    Return user-specified model path, or the newest *with_HMPV_VBOF.xml under OUTPUT_DIR.
    """
    if explicit is not None:
        return explicit
    matches = sorted(
        OUTPUT_DIR.glob(f"*{INTEGRATED_MODEL_XML_SUFFIX}"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError(
            f"No integrated model (*{INTEGRATED_MODEL_XML_SUFFIX}) under {OUTPUT_DIR}. "
            "Run integrate_model.py or pass --model."
        )
    chosen = matches[0]
    logger.info("Using integrated model: %s", chosen)
    return chosen


def load_target_tables(
    reaction_csv: Path,
    gene_csv: Path,
) -> pd.DataFrame:
    """
    Load selective reaction and gene CSVs and return one combined table with:
    - target_id: reaction_id or gene_id
    - knockout_type: 'reaction' | 'gene'
    """
    frames: List[pd.DataFrame] = []

    if reaction_csv.exists():
        rdf = pd.read_csv(reaction_csv)
        if not rdf.empty and "reaction_id" in rdf.columns:
            frames.append(
                pd.DataFrame(
                    {
                        "target_id": rdf["reaction_id"].astype(str),
                        "knockout_type": "reaction",
                    }
                )
            )
        elif not rdf.empty:
            logger.warning("%s has no reaction_id column; skipping.", reaction_csv)
    else:
        logger.warning("Reaction targets file not found: %s", reaction_csv)

    if gene_csv.exists():
        gdf = pd.read_csv(gene_csv)
        if not gdf.empty and "gene_id" in gdf.columns:
            frames.append(
                pd.DataFrame(
                    {
                        "target_id": gdf["gene_id"].astype(str),
                        "knockout_type": "gene",
                    }
                )
            )
        elif not gdf.empty:
            logger.warning("%s has no gene_id column; skipping.", gene_csv)
    else:
        logger.warning("Gene targets file not found: %s", gene_csv)

    if not frames:
        return pd.DataFrame(
            columns=["target_id", "knockout_type"]
        )

    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["target_id", "knockout_type"]).reset_index(drop=True)
    return out


def _resolve_biomass_reaction(model: Model, reaction_id: str):
    """Return reaction if present; try common SBML prefix variant."""
    if reaction_id in model.reactions:
        return model.reactions.get_by_id(reaction_id)
    alt = f"R_{reaction_id}"
    if alt in model.reactions:
        return model.reactions.get_by_id(alt)
    raise KeyError(f"Biomass reaction {reaction_id!r} (or {alt!r}) not in model.")


def _fva_min_max(
    model: Model,
    reaction_id: str,
    fraction: float,
) -> Tuple[float, float]:
    """
    Run FVA for a single reaction with the *current* model objective.

    COBRApy maximizes the objective, then constrains it to >= fraction * optimum
    and computes minimum/maximum flux for each reaction in reaction_list.
    """
    rxn = model.reactions.get_by_id(reaction_id)
    fva = flux_variability_analysis(
        model,
        reaction_list=[rxn],
        fraction_of_optimum=fraction,
    )
    # Index is reaction id; columns are 'minimum', 'maximum'
    row = fva.loc[reaction_id]
    return float(row["minimum"]), float(row["maximum"])


def run_knockout_fva_row(
    model: Model,
    knockout_type: str,
    target_id: str,
    vbof_id: str,
    host_bof_id: str,
    fraction: float,
    run_host: bool,
) -> dict:
    """
    Apply one knockout (reaction or gene), then:
    1) Maximize VBOF; FVA on VBOF at fraction_of_optimum.
    2) Optionally maximize host BOF; FVA on host BOF at fraction_of_optimum.

    All changes occur inside ``with model:`` (caller must open the context).
    """
    # --- Apply knockout -------------------------------------------------
    if knockout_type == "reaction":
        if target_id not in model.reactions:
            raise KeyError(f"Reaction {target_id!r} not in model.")
        model.reactions.get_by_id(target_id).knock_out()
    elif knockout_type == "gene":
        if target_id not in model.genes:
            raise KeyError(f"Gene {target_id!r} not in model.")
        model.genes.get_by_id(target_id).knock_out()
    else:
        raise ValueError(f"Unknown knockout_type: {knockout_type!r}")

    # --- Viral biomass: set objective, FVA on HMPV_VBOF -----------------
    vbof_rxn = model.reactions.get_by_id(vbof_id)
    model.objective = vbof_rxn
    virus_min, virus_max = _fva_min_max(model, vbof_id, fraction)

    host_min = host_max = float("nan")
    if run_host:
        # --- Host biomass: switch objective, FVA on host BOF ----------
        host_rxn = _resolve_biomass_reaction(model, host_bof_id)
        model.objective = host_rxn
        host_min, host_max = _fva_min_max(model, host_bof_id, fraction)

    return {
        "target_id": target_id,
        "virus_fva_min": virus_min,
        "virus_fva_max": virus_max,
        "host_fva_min": host_min,
        "host_fva_max": host_max,
    }


def analyze_targets(
    model: Model,
    targets: pd.DataFrame,
    vbof_id: str,
    host_bof_id: str,
    fraction: float,
    run_host: bool,
) -> pd.DataFrame:
    """Loop targets; each uses ``with model:`` so the base model resets."""
    rows = []
    for i, row in targets.iterrows():
        tid = row["target_id"]
        ktype = row["knockout_type"]
        logger.info("(%s/%s) %s knockout: %s", i + 1, len(targets), ktype, tid)
        try:
            # COBRApy History: entering reverts on exit, so knockouts never accumulate
            with model:
                rec = run_knockout_fva_row(
                    model,
                    knockout_type=ktype,
                    target_id=tid,
                    vbof_id=vbof_id,
                    host_bof_id=host_bof_id,
                    fraction=fraction,
                    run_host=run_host,
                )
            rec["knockout_type"] = ktype
            rec["status"] = "ok"
            rows.append(rec)
        except Exception as exc:  # noqa: BLE001 — log and continue per target
            logger.exception("Failed for %s %s: %s", ktype, tid, exc)
            rows.append(
                {
                    "target_id": tid,
                    "knockout_type": ktype,
                    "virus_fva_min": np.nan,
                    "virus_fva_max": np.nan,
                    "host_fva_min": np.nan,
                    "host_fva_max": np.nan,
                    "status": f"error: {exc}",
                }
            )

    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--model",
        type=Path,
        default=None,
        help="Integrated SBML model path (default: newest under output_HBEC)",
    )
    p.add_argument(
        "--reaction-csv",
        type=Path,
        default=SELECTIVE_REACTION_TARGETS_PATH,
        help="selective_reaction_targets.csv",
    )
    p.add_argument(
        "--gene-csv",
        type=Path,
        default=SELECTIVE_TARGETS_PATH,
        help="selective_antiviral_targets.csv",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=SELECTIVE_TARGETS_FVA_PATH,
        help="Output CSV path (default: dual_objective_analysis/selective_targets_fva.csv)",
    )
    p.add_argument(
        "--fraction",
        type=float,
        default=0.9,
        help="FVA fraction_of_optimum (default: 0.9)",
    )
    p.add_argument(
        "--no-host-fva",
        action="store_true",
        help="Only run FVA with VBOF objective; skip host BOF FVA",
    )
    p.add_argument("--vbof-id", type=str, default=VBOF_REACTION_ID)
    p.add_argument("--host-bof-id", type=str, default=HOST_BOF_REACTION_ID)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    logger.info("=" * 70)
    logger.info("SELECTIVE TARGETS FLUX VARIABILITY ANALYSIS")
    logger.info("=" * 70)

    # =========================================================================
    # STEP 1: Locate integrated host + HMPV_VBOF SBML and load once
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 1: LOADING INTEGRATED MODEL")
    logger.info("=" * 70)

    model_path = resolve_integrated_model_path(args.model)
    model = load_integrated_model(model_path, vbof_id=args.vbof_id)

    # =========================================================================
    # STEP 2: Verify biomass reactions exist before the knockout loop
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: VALIDATING BIOMASS REACTIONS")
    logger.info("=" * 70)

    if args.vbof_id not in model.reactions:
        raise KeyError(f"VBOF reaction {args.vbof_id!r} not in model.")
    _resolve_biomass_reaction(model, args.host_bof_id)
    logger.info("VBOF reaction: %s — OK", args.vbof_id)
    logger.info("Host BOF reaction: %s — OK", args.host_bof_id)

    # =========================================================================
    # STEP 3: Build target queue from reaction- and gene-level CSVs
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: LOADING TARGET TABLES")
    logger.info("=" * 70)

    targets = load_target_tables(args.reaction_csv, args.gene_csv)
    if targets.empty:
        logger.error("No targets loaded; exiting.")
        sys.exit(1)
    logger.info(
        "Loaded %s targets (%s reaction, %s gene)",
        len(targets),
        (targets["knockout_type"] == "reaction").sum(),
        (targets["knockout_type"] == "gene").sum(),
    )

    # =========================================================================
    # STEP 4: Run FVA for each target (``with model:`` resets after each)
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 4: RUNNING KNOCKOUT FVA")
    logger.info("=" * 70)
    logger.info(
        "fraction_of_optimum=%.2f | host FVA: %s",
        args.fraction,
        "enabled" if not args.no_host_fva else "disabled",
    )

    result = analyze_targets(
        model,
        targets,
        vbof_id=args.vbof_id,
        host_bof_id=args.host_bof_id,
        fraction=args.fraction,
        run_host=not args.no_host_fva,
    )

    # =========================================================================
    # STEP 5: Write output CSV
    # =========================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 5: SAVING RESULTS")
    logger.info("=" * 70)

    col_order = [
        "target_id",
        "knockout_type",
        "virus_fva_min",
        "virus_fva_max",
        "host_fva_min",
        "host_fva_max",
        "status",
    ]
    result = result[[c for c in col_order if c in result.columns]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    logger.info("Wrote %s (%s rows)", args.output, len(result))

    ok_count = (result["status"] == "ok").sum()
    err_count = len(result) - ok_count
    logger.info("Successful: %s | Errors: %s", ok_count, err_count)

    logger.info("\n" + "=" * 70)
    logger.info("SELECTIVE TARGETS FVA COMPLETE")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
