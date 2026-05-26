"""Flag systematic Visium outlier spots by barcode matching."""

import numpy as np
import pandas as pd
import os


def _load_biased_spots() -> pd.DataFrame:
    """Load the built-in biased_spots dataset."""
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    csv_path = os.path.join(data_dir, "biased_spots.csv")
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)

    # Fallback: load from the R-generated fixture
    fixture_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    csv_path = os.path.join(fixture_dir, "biased_spots.csv")
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)

    raise FileNotFoundError(
        "biased_spots.csv not found. Ensure the data file is in the package."
    )


def flag_visium_outliers(
    adata,
) -> "anndata.AnnData":
    """Flag systematic Visium outlier spots based on array coordinates.

    Matches spot array coordinates against a known list of biased spots
    (technical outliers in Visium arrays). Adds a boolean
    ``systematic_outliers`` column to ``adata.obs``.

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data with ``array_row`` and ``array_col``
        columns in ``adata.obs``.

    Returns
    -------
    anndata.AnnData
        Modified AnnData with ``systematic_outliers`` column in ``adata.obs``.
    """
    if "array_row" not in adata.obs.columns:
        raise ValueError("'array_row' not found in adata.obs.")
    if "array_col" not in adata.obs.columns:
        raise ValueError("'array_col' not found in adata.obs.")

    biased_spots = _load_biased_spots()

    drop_mask = np.zeros(len(adata), dtype=bool)

    for _, row in biased_spots.iterrows():
        row_match = adata.obs["array_row"].values == row["row"]
        col_match = adata.obs["array_col"].values == row["col"]
        drop_mask = drop_mask | (row_match & col_match)

    adata.obs["systematic_outliers"] = drop_mask

    return adata
