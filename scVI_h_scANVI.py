"""
Train scVI on raw counts (layer 'counts'), then scANVI on cell_type and on h_cell_type.
Requires GPU for practical runtime; CPU is supported but slow.
"""

from __future__ import annotations

import sys

import scanpy as sc
import scvi
import torch
from rich import print

import pipeline_paths as pp

scvi.settings.seed = 0
torch.manual_seed(0)

INPUT_FILE = pp.path_combined_pre_integration()
OUTPUT_FILE_1 = pp.path_combined_post_scvi()
OUTPUT_FILE_2 = pp.path_combined_post_scanvi()
SCANVI_MODEL_SAVE = pp.path_scanvi_model()

if not INPUT_FILE.is_file():
    print(f"ERROR: Input file not found: {INPUT_FILE}")
    sys.exit(1)

if torch.cuda.is_available():
    print(f"  GPU available: {torch.cuda.get_device_name(0)}")
    print(f"  CUDA version: {torch.version.cuda}")
else:
    print("  WARNING: No GPU detected, training will use CPU (may be slow)")

print(f"  Loading combined dataset from: {INPUT_FILE}")
try:
    adata = sc.read_h5ad(INPUT_FILE)
    print(f"  Loaded dataset with shape: {adata.shape}")
except Exception as e:
    print(f"ERROR: Failed to load data: {e}")
    sys.exit(1)

print("Training scVI and scANVI on cell_type labels")

print("  Initializing scVI model")
try:
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor_id_tissue")
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    print("  Model initialized successfully")
except Exception as e:
    print(f"ERROR: Failed to initialize model: {e}")
    sys.exit(1)

print("  Training scVI model")
try:
    model.train(max_epochs=400, check_val_every_n_epoch=50)
    print("  scVI training completed")
except Exception as e:
    print(f"ERROR: Training failed: {e}")
    sys.exit(1)

print("  Saving scVI latent space")
try:
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    OUTPUT_FILE_1.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUTPUT_FILE_1)
    print("  scVI latent space saved to X_scVI latent key")
except Exception as e:
    print(f"ERROR: Failed to save scVI latent space: {e}")
    sys.exit(1)

print("  Initializing scANVI model")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown",
)
print("  scANVI model initialized successfully")
print("  Training scANVI model")
try:
    scanvi_model.train(max_epochs=20, n_samples_per_label=100)
    print("  scANVI training completed")
except Exception as e:
    print(f"ERROR: Training failed: {e}")
    sys.exit(1)

print("  Saving scANVI latent space")
try:
    SCANVI_LATENT_KEY = "X_scANVI"
    adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
    print("  scANVI latent space saved to X_scANVI latent key")
except Exception as e:
    print(f"ERROR: Failed to save scANVI latent space: {e}")
    sys.exit(1)

print("  scANVI pipeline round 1/2 completed successfully!")

print("Training scANVI on h_cell_type labels")

print("  Initializing scANVI model")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="h_cell_type",
    unlabeled_category="Unknown",
)
print("  scANVI model initialized successfully")
print("  Training scANVI model")
try:
    scanvi_model.train(max_epochs=20, n_samples_per_label=100)
    print("  scANVI training completed")
except Exception as e:
    print(f"ERROR: Training failed: {e}")
    sys.exit(1)

print(f"  Saving SCANVI model to: {SCANVI_MODEL_SAVE}")
try:
    SCANVI_MODEL_SAVE.parent.mkdir(parents=True, exist_ok=True)
    scanvi_model.save(str(SCANVI_MODEL_SAVE), overwrite=True)
    print("  SCANVI model saved successfully")
except Exception as e:
    print(f"ERROR: Failed to save SCANVI model: {e}")
    sys.exit(1)

print("  Saving h_scANVI latent space")
try:
    H_SCANVI_LATENT_KEY = "X_h_scANVI"
    adata.obsm[H_SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
    adata.write_h5ad(OUTPUT_FILE_2)
    print(
        f"  scANVI latent space saved to X_h_scANVI latent key and written to: {OUTPUT_FILE_2}"
    )
except Exception as e:
    print(f"ERROR: Failed to save scANVI latent space and write to file: {e}")
    sys.exit(1)

print("  scANVI pipeline round 2/2 completed successfully!")
