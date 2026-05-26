"""Shared fixtures for SpotSweeper parity tests."""

import json
import os
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def load_visium_fixture():
    """Load the Visium DLPFC fixture as AnnData."""
    import anndata as ad

    counts = mmread(os.path.join(DATA_DIR, "fixture_counts.mtx"))
    if hasattr(counts, "toarray"):
        counts = counts.T  # transpose to spots x genes
    else:
        counts = counts.T

    coords = pd.read_csv(
        os.path.join(DATA_DIR, "fixture_spatial_coords.csv"), index_col=0
    )
    metadata = pd.read_csv(
        os.path.join(DATA_DIR, "fixture_metadata.csv"), index_col=0
    )
    features = pd.read_csv(os.path.join(DATA_DIR, "feature_names.csv"))
    spots = pd.read_csv(os.path.join(DATA_DIR, "spot_names.csv"))

    adata = ad.AnnData(
        X=csr_matrix(counts) if not hasattr(counts, "toarray") else csr_matrix(counts),
        obs=metadata,
        var=pd.DataFrame(index=features.iloc[:, 0].values),
    )
    adata.obs_names = spots.iloc[:, 0].values
    adata.obsm["spatial"] = coords.values

    return adata


def load_artifact_fixture():
    """Load the DLPFC artifact fixture as AnnData."""
    import anndata as ad

    counts = mmread(os.path.join(DATA_DIR, "fixture_artifact_counts.mtx"))
    if hasattr(counts, "toarray"):
        counts = counts.T
    else:
        counts = counts.T

    coords = pd.read_csv(
        os.path.join(DATA_DIR, "fixture_artifact_coords.csv"), index_col=0
    )
    metadata = pd.read_csv(
        os.path.join(DATA_DIR, "fixture_artifact_metadata.csv"), index_col=0
    )
    features = pd.read_csv(os.path.join(DATA_DIR, "fixture_artifact_features.csv"))
    spots = pd.read_csv(os.path.join(DATA_DIR, "fixture_artifact_spots.csv"))

    adata = ad.AnnData(
        X=csr_matrix(counts),
        obs=metadata,
        var=pd.DataFrame(index=features.iloc[:, 0].values),
    )
    adata.obs_names = spots.iloc[:, 0].values
    adata.obsm["spatial"] = coords.values

    return adata


def load_r_knn_indices():
    """Load R kNN indices (0-based)."""
    r_knn = pd.read_csv(os.path.join(DATA_DIR, "r_knn_indices_k36.csv")).values
    return r_knn - 1  # Convert from 1-based to 0-based


def load_reference_outputs():
    """Load R reference outputs."""
    with open(os.path.join(DATA_DIR, "reference_outputs.json")) as f:
        return json.load(f)


def load_artifact_reference_outputs():
    """Load R artifact reference outputs."""
    with open(os.path.join(DATA_DIR, "reference_artifact_outputs.json")) as f:
        return json.load(f)
