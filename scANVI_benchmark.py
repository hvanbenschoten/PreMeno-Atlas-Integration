"""
Integration benchmarking with scib-metrics (GPU k-NN via cuML optional).

Expects the downsampled object from notebooks/downsample_500k.ipynb with latent keys
X_scVI, X_scANVI, X_h_scANVI (or compatible obsm) and X_pca for unintegrated baseline.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import scanpy as sc
from rich import print
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from scib_metrics.nearest_neighbors import NeighborsResults

import pipeline_paths as pp

INPUT_FILE = pp.path_combined_500k()
OUTPUT_DIR = pp.INTEGRATION_OUTPUT

if not INPUT_FILE.is_file():
    print(f"ERROR: Input file not found: {INPUT_FILE}")
    print("Run notebooks/downsample_500k.ipynb first (or set ATLAS_* paths).")
    sys.exit(1)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output directory: {OUTPUT_DIR}")

try:
    adata = sc.read_h5ad(INPUT_FILE)
    print(f"  Loaded dataset with shape: {adata.shape}")
except Exception as e:
    print(f"ERROR: Failed to load data: {e}")
    sys.exit(1)

print("  Preprocessing data for benchmark...")
adata.obsm["scVI"] = adata.obsm["X_scVI"]
adata.obsm["scANVI"] = adata.obsm["X_scANVI"]
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]


def cuml_nearest_neighbors(X: np.ndarray, k: int):
    """GPU k-NN via RAPIDS cuML (optional)."""
    from cuml.neighbors import NearestNeighbors

    X = np.ascontiguousarray(X, dtype=np.float32)
    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    indices = np.asarray(indices, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.float64)
    return NeighborsResults(indices=indices, distances=distances)


def _get_neighbor_computer():
    try:
        from cuml.neighbors import NearestNeighbors  # noqa: F401

        return cuml_nearest_neighbors
    except ImportError:
        print("  WARNING: cuML not found. Using default CPU nearest neighbors.")
        return None


biocons = BioConservation(isolated_labels=False)

print("  Benchmarking...")
start = time.time()
bm = Benchmarker(
    adata,
    batch_key="donor_id_tissue",
    label_key="cell_type",
    embedding_obsm_keys=["Unintegrated", "scANVI"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,
)
neighbor_computer = _get_neighbor_computer()
if neighbor_computer is not None:
    bm.prepare(neighbor_computer=neighbor_computer)
else:
    bm.prepare()
bm.benchmark()

print("  Plotting benchmark results...")
tab_scaled = bm.plot_results_table(show=False)
scaled_pdf_path = OUTPUT_DIR / "scANVI_benchmark_results.pdf"
tab_scaled.figure.savefig(scaled_pdf_path, bbox_inches="tight", dpi=300)
print(f"  ✓ Saved scaled results to: {scaled_pdf_path}")

tab_unscaled = bm.plot_results_table(min_max_scale=False, show=False)
unscaled_pdf_path = OUTPUT_DIR / "scANVI_benchmark_results_unscaled.pdf"
tab_unscaled.figure.savefig(unscaled_pdf_path, bbox_inches="tight", dpi=300)
print(f"  ✓ Saved unscaled results to: {unscaled_pdf_path}")

df = bm.get_results(min_max_scale=False)
csv_path = OUTPUT_DIR / "scANVI_benchmark_results_unscaled.csv"
df.to_csv(csv_path, index=False)
print(f"  ✓ Saved results table to: {csv_path}")

print("\n  Benchmarking completed successfully!")
print(f"  Total time: {time.time() - start:.2f} seconds")
