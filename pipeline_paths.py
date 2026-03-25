"""
Central path configuration for the integration pipeline.

Override defaults with environment variables (see README) or edit paths here.
Repository root = directory containing this file (the clone root).
"""

from __future__ import annotations

import os
from pathlib import Path


def repo_root() -> Path:
    env = os.environ.get("ATLAS_PIPELINE_ROOT", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parent


ROOT = repo_root()
DATA = Path(os.environ.get("ATLAS_DATA_DIR", ROOT / "data")).expanduser().resolve()

# Per-stage directories (created on demand by scripts)
PREPROCESSED_H5AD = Path(
    os.environ.get("ATLAS_PREPROCESSED_H5AD", DATA / "preprocessed_h5ad")
).expanduser().resolve()
PROCESSED_H5AD = Path(
    os.environ.get("ATLAS_PROCESSED_H5AD", DATA / "processed_h5ad")
).expanduser().resolve()
INTEGRATION_OUTPUT = Path(
    os.environ.get("ATLAS_INTEGRATION_OUTPUT", DATA / "integration_output")
).expanduser().resolve()

# Harmonization table: columns tissue_id, cell_type, h_cell_type
HARMONIZATION_KEY_CSV = Path(
    os.environ.get("ATLAS_HARMONIZATION_KEY", ROOT / "config" / "harmonization_key.csv")
).expanduser().resolve()

# Optional: one tissue name per line (see config/tissues.txt)
TISSUES_FILE = Path(
    os.environ.get("ATLAS_TISSUES_FILE", ROOT / "config" / "tissues.txt")
).expanduser().resolve()

# --- Intermediate h5ad filenames (under INTEGRATION_OUTPUT) ---
COMBINED_PRE_INTEGRATION = "combined_pre_integration.h5ad"
COMBINED_POST_SCVI = "combined_pre_integration_scVI.h5ad"
COMBINED_POST_SCANVI = "combined_integrated_h_scANVI.h5ad"
COMBINED_ANALYZED = "combined_integrated_h_scANVI_analyzed.h5ad"
COMBINED_500K = "combined_integrated_h_scANVI_analyzed_500k.h5ad"
SCANVI_MODEL_DIR_NAME = "scanvi_h_cell_type_model"


def path_combined_pre_integration() -> Path:
    return INTEGRATION_OUTPUT / COMBINED_PRE_INTEGRATION


def path_combined_post_scvi() -> Path:
    return INTEGRATION_OUTPUT / COMBINED_POST_SCVI


def path_combined_post_scanvi() -> Path:
    return INTEGRATION_OUTPUT / COMBINED_POST_SCANVI


def path_combined_analyzed() -> Path:
    return INTEGRATION_OUTPUT / COMBINED_ANALYZED


def path_combined_500k() -> Path:
    return INTEGRATION_OUTPUT / COMBINED_500K


def path_scanvi_model() -> Path:
    return INTEGRATION_OUTPUT / SCANVI_MODEL_DIR_NAME
