# Multi-tissue scRNA-seq integration pipeline (scVI / scANVI)

Reference workflow for harmonizing CellxGene-style and custom `h5ad` objects, building a concatenated atlas with harmonized cell-type labels, training **scVI** and **scANVI**, exploring results, downsampling for benchmarking, and running **scib-metrics** benchmarks.

This repository is **not** a pip-installable package: run notebooks interactively and execute Python scripts from the repository root (or set paths as below).

## Pipeline overview

| Step | Artifact | Description |
|------|----------|-------------|
| 1 | `notebooks/preprocessing.ipynb` | Per-dataset QC and processing from raw / CellxGene-style `h5ad` into `data/preprocessed_h5ad/` → `data/processed_h5ad/`. |
| 2 | `notebooks/h5ad_inspect.ipynb` | Compare object structure (e.g. curated vs CellxGene) and fix inconsistencies. |
| 3 | `preintegration.py` | Concatenate tissues (using `.raw` counts), merge `h_cell_type` from a CSV key, HVG selection, write `data/integration_output/combined_pre_integration.h5ad`. |
| 4 | `scVI_h_scANVI.py` | Train scVI + two rounds of scANVI (`cell_type`, then `h_cell_type`); save latents and final SCANVI model. |
| 5 | `notebooks/h_scANVI_results.ipynb` | Load integrated object (`X_h_scANVI` focus), QC latent spaces, UMAPs, optional analyses. |
| 6 | `notebooks/downsample_500k.ipynb` | Memory-safe downsample (backed mode) for metrics. |
| 7 | `scANVI_benchmark.py` | scib-metrics benchmark (GPU k-NN optional via RAPIDS cuML). |

## Layout

```text
integration_pipeline/
  pipeline_paths.py       # Central paths (edit or use env vars)
  preintegration.py
  scVI_h_scANVI.py
  scANVI_benchmark.py
  requirements.txt
  README.md
  config/
    tissues.txt                    # Stems of `<stem>.h5ad` in processed_h5ad
    harmonization_key.example.csv  # Copy to harmonization_key.csv and complete
  data/                            # Created locally; h5ad files are gitignored
    preprocessed_h5ad/
    processed_h5ad/
    integration_output/
  notebooks/
    preprocessing.ipynb
    h5ad_inspect.ipynb
    h_scANVI_results.ipynb
    downsample_500k.ipynb
  slurm/
    run_preintegration.sbatch
    run_scANVI.sbatch
    run_scANVI_benchmark.sbatch
```

## Setup

1. Clone this directory as its own Git repository (or copy the folder into your project).

2. Create a conda (or venv) environment and install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

   Optional notebook pieces: **scrublet** (doublets), **rpy2** + R (only if you use the R-backed cells in `preprocessing.ipynb`). For **GPU k-NN** in `scANVI_benchmark.py`, install **RAPIDS cuML** matching your CUDA stack (conda is often easiest).

3. **Paths:** By default, all data live under `./data/...` relative to the repository root (the directory that contains `pipeline_paths.py`). To point elsewhere:

   | Variable | Meaning |
   |----------|---------|
   | `ATLAS_PIPELINE_ROOT` | Repository root (folder with `pipeline_paths.py`) |
   | `ATLAS_DATA_DIR` | Parent of `preprocessed_h5ad`, `processed_h5ad`, `integration_output` (default: `<root>/data`) |
   | `ATLAS_PREPROCESSED_H5AD` | Input folder for step 1 outputs |
   | `ATLAS_PROCESSED_H5AD` | Per-tissue processed `h5ad` for step 3 |
   | `ATLAS_INTEGRATION_OUTPUT` | Combined objects and benchmarks |
   | `ATLAS_HARMONIZATION_KEY` | CSV path (default: `config/harmonization_key.csv`) |
   | `ATLAS_TISSUES_FILE` | List of tissue stems (default: `config/tissues.txt`) |

4. **Harmonization key:** Copy `config/harmonization_key.example.csv` to `config/harmonization_key.csv`. Required columns: `tissue_id`, `cell_type`, `h_cell_type`. Rows should cover the `(cell_type, tissue_id)` pairs present in your combined object.

5. **Jupyter:** Start Jupyter from the repository root **or** set `ATLAS_PIPELINE_ROOT` before launching. Each notebook begins with a short cell that imports `pipeline_paths` after locating the repo.

6. **SLURM:** Edit `ATLAS_PIPELINE_ROOT`, partition names, memory, and `conda activate` lines in `slurm/*.sbatch` for your cluster.

## File naming (integration output)

Scripts and notebooks assume these filenames under `integration_output` (configurable in `pipeline_paths.py`):

- `combined_pre_integration.h5ad` — after `preintegration.py`
- `combined_pre_integration_scVI.h5ad` — after scVI latent write
- `combined_integrated_h_scANVI.h5ad` — after full `scVI_h_scANVI.py`
- `combined_integrated_h_scANVI_analyzed.h5ad` — after `h_scANVI_results.ipynb`
- `combined_integrated_h_scANVI_analyzed_500k.h5ad` — after `downsample_500k.ipynb`

## Notes

- **HVG selection** in `preintegration.py` uses `flavor="seurat"` with `batch_key="donor_id_tissue"` to avoid LOESS issues on very large batch counts; adjust if your atlas is smaller.
- **scVI** expects raw counts in `adata.layers["counts"]` (written during preintegration).
- **Benchmark script** maps `X_scVI` / `X_scANVI` into `obsm` keys expected by the benchmarker; ensure your downsampled object still contains `X_pca` for the unintegrated baseline.

## License

Add a `LICENSE` file appropriate to your lab or publication policy before publishing on GitHub.
