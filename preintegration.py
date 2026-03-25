#!/usr/bin/env python3
"""
Pre-integration: load tissue h5ad objects, harmonize cell types, HVG selection, save combined h5ad.

Linear counterpart to exploratory pre_integration notebooks. Run after per-tissue processing
(`notebooks/preprocessing.ipynb`, `notebooks/h5ad_inspect.ipynb`).
"""

from __future__ import annotations

import sys
from typing import List

import pandas as pd
import scanpy as sc

import pipeline_paths as pp


def _load_tissue_list() -> List[str]:
    path = pp.TISSUES_FILE
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing tissues list: {path}\n"
            "Copy config/tissues.txt from the example or create your own (one stem per line)."
        )
    lines: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        lines.append(line)
    if not lines:
        raise ValueError(f"No tissue names found in {path}")
    return lines


PROCESSED_TISSUES: List[str] = []  # filled in main() from tissues file


def load_full_tissues() -> sc.AnnData:
    """
    Load each processed tissue h5ad from PROCESSED_H5AD, use .raw for counts,
    standardize metadata, concatenate.
    """
    print("=== Load and combine processed tissues (using .raw = raw counts) ===")
    adatas = []
    missing_raw = []
    input_dir = pp.PROCESSED_H5AD

    for tissue in PROCESSED_TISSUES:
        path = input_dir / f"{tissue}.h5ad"
        if not path.is_file():
            print(f"  ⚠ File not found, skipping: {path}")
            continue

        print(f"Loading {tissue} from {path} ...")
        adata = sc.read_h5ad(path)
        print(f"  Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

        if adata.raw is None:
            missing_raw.append(tissue)
            continue

        adata_raw = sc.AnnData(
            X=adata.raw.X.copy(),
            var=adata.raw.var.copy(),
            obs=adata.obs.copy(),
        )

        adata_raw.obs["tissue"] = tissue
        adata_raw.obs["tissue_id"] = tissue

        donor_col = None
        for col in adata_raw.obs.columns:
            if "donor" in str(col).lower():
                donor_col = col
                break

        if donor_col is not None:
            adata_raw.obs["donor_id"] = adata_raw.obs[donor_col].astype(str)
        else:
            print(f"  ⚠ No donor column found in {tissue}; using '{tissue}_unknown'")
            adata_raw.obs["donor_id"] = f"{tissue}_unknown"

        adata_raw.obs["donor_id_tissue"] = (
            adata_raw.obs["donor_id"].astype(str)
            + "_"
            + adata_raw.obs["tissue_id"].astype(str)
        )

        adatas.append(adata_raw)

    if missing_raw:
        raise RuntimeError(
            f"Tissues missing .raw attribute (required for raw-count concatenation): {missing_raw}. "
            "Ensure all tissue h5ad files have .raw with raw counts."
        )

    if not adatas:
        raise RuntimeError("No tissue files could be loaded; nothing to combine.")

    print("\nConcatenating all tissues from .raw (sc.concat) ...")
    combined = sc.concat(
        adatas,
        label="tissue",
        index_unique="-",
        join="outer",
    )

    print("\nCombined dataset summary:")
    print(f"  Total cells: {combined.n_obs:,}")
    print(f"  Total genes: {combined.n_vars:,}")
    print(f"  Tissues: {combined.obs['tissue'].nunique()}")
    if "donor_id" in combined.obs.columns:
        print(f"  Donors: {combined.obs['donor_id'].nunique()}")
    if "donor_id_tissue" in combined.obs.columns:
        print(
            "  Unique donor-tissue combinations: "
            f"{combined.obs['donor_id_tissue'].nunique():,}"
        )
    if "cell_type" in combined.obs.columns:
        print(f"  Cell types: {combined.obs['cell_type'].nunique()}")

    return combined


def add_harmonized_cell_types(combined: sc.AnnData) -> sc.AnnData:
    """Merge harmonization key; drop cells with unknown harmonized labels."""
    key_path = pp.HARMONIZATION_KEY_CSV
    print("\n=== Add harmonized cell type labels (h_cell_type) ===")
    if not key_path.is_file():
        raise FileNotFoundError(
            f"Harmonization key CSV not found: {key_path}\n"
            "Copy config/harmonization_key.example.csv to config/harmonization_key.csv and edit, "
            "or set ATLAS_HARMONIZATION_KEY."
        )

    harmonization_key = pd.read_csv(key_path)
    required_cols = {"cell_type", "tissue_id", "h_cell_type"}
    missing = required_cols - set(harmonization_key.columns)
    if missing:
        raise ValueError(
            f"Harmonization key missing required columns: {sorted(missing)}"
        )

    print(
        f"  Loaded harmonization key with "
        f"{len(harmonization_key):,} rows from:\n    {key_path}"
    )
    print(
        "  Unique (cell_type, tissue_id) pairs in key: "
        f"{harmonization_key[['cell_type', 'tissue_id']].drop_duplicates().shape[0]:,}"
    )

    if "cell_type" not in combined.obs.columns or "tissue_id" not in combined.obs.columns:
        raise ValueError(
            "combined.obs must contain 'cell_type' and 'tissue_id' before harmonization."
        )

    combined_obs_df = combined.obs[["cell_type", "tissue_id"]].copy()
    combined_obs_df = combined_obs_df.reset_index().rename(columns={"index": "obs_index"})

    key_sub = harmonization_key[["cell_type", "tissue_id", "h_cell_type"]].copy()
    merged = combined_obs_df.merge(
        key_sub,
        on=["cell_type", "tissue_id"],
        how="left",
    )

    combined.obs["h_cell_type"] = merged["h_cell_type"].values

    n_total = combined.n_obs
    n_missing = combined.obs["h_cell_type"].isna().sum()
    print(
        f"  Cells with missing harmonized label (h_cell_type is NaN): "
        f"{n_missing:,} / {n_total:,}"
    )

    print("\n=== Remove cells with unknown / missing harmonized labels ===")
    unknown_mask = combined.obs["h_cell_type"].isna() | (
        combined.obs["h_cell_type"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("unknown")
    )

    n_unknown = int(unknown_mask.sum())
    print(f"  Cells with 'Unknown'/missing h_cell_type: {n_unknown:,}")

    if n_unknown > 0:
        before = combined.n_obs
        combined = combined[~unknown_mask].copy()
        after = combined.n_obs
        print(
            f"  Removed {n_unknown:,} unknown cells; "
            f"{before:,} -> {after:,} cells"
        )
    else:
        print("  No unknown cells found; nothing removed.")

    print(
        f"  Unique harmonized cell types (h_cell_type): "
        f"{combined.obs['h_cell_type'].nunique():,}"
    )
    if "cell_type" in combined.obs.columns:
        print(
            f"  Unique original cell types (cell_type): "
            f"{combined.obs['cell_type'].nunique():,}"
        )

    return combined


def perform_hvg_selection(combined: sc.AnnData) -> sc.AnnData:
    """Normalize+log for HVG detection; subset to HVGs; keep counts in layers['counts']."""
    print("\n=== Highly variable gene selection ===")

    if "donor_id_tissue" not in combined.obs.columns:
        raise ValueError(
            "combined.obs must contain 'donor_id_tissue' for batch-aware HVG selection."
        )

    print("  Saving raw counts to combined.layers['counts'] ...")
    combined.layers["counts"] = combined.X.copy()

    print("  Normalizing and log-transforming X for HVG detection ...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    print("  Computing highly variable genes (HVGs) ...")
    sc.pp.highly_variable_genes(
        combined,
        n_top_genes=4000,
        batch_key="donor_id_tissue",
        flavor="seurat",
    )

    n_hvg = int(combined.var.get("highly_variable", pd.Series(False)).sum())
    print(f"  Selected {n_hvg:,} HVGs")

    if n_hvg == 0:
        raise RuntimeError("No highly variable genes were selected.")

    combined = combined[:, combined.var["highly_variable"]].copy()
    print(
        f"  After HVG subset: {combined.n_obs:,} cells, {combined.n_vars:,} genes"
    )

    combined.raw = combined

    return combined


def main() -> None:
    global PROCESSED_TISSUES
    PROCESSED_TISSUES = _load_tissue_list()

    pp.INTEGRATION_OUTPUT.mkdir(parents=True, exist_ok=True)
    output_h5ad = pp.path_combined_pre_integration()

    combined = load_full_tissues()
    combined = add_harmonized_cell_types(combined)
    combined = perform_hvg_selection(combined)

    print("\n=== Save integration-ready object ===")
    for col in combined.obs.columns:
        if combined.obs[col].dtype == object or str(combined.obs[col].dtype) == "category":
            combined.obs[col] = combined.obs[col].astype(str)
    print(f"  Writing: {output_h5ad}")
    combined.write_h5ad(output_h5ad)
    print(
        f"  ✓ Saved: {output_h5ad}\n"
        f"  Final: {combined.n_obs:,} cells, {combined.n_vars:,} genes"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n[ERROR] preintegration.py failed: {e}", file=sys.stderr)
        raise
