"""Local outlier detection using kNN-based modified z-scores."""

from typing import Optional
import numpy as np
import pandas as pd
from spotsweeper.utils import find_knn



def _modified_zscore(x: np.ndarray, s: float = 1.4826) -> np.ndarray:
    """Compute modified z-scores using MAD (matching spatialEco::outliers).

    Formula: 0.6745 * (x - median(x)) / mad(x)
    where mad(x) = s * median(|x - median(x)|)

    Parameters
    ----------
    x : np.ndarray
        Input values.
    s : float
        Scaling constant for MAD (default 1.4826, consistent for normal data).

    Returns
    -------
    np.ndarray
        Modified z-scores.
    """
    median_x = np.median(x)
    mad = s * np.median(np.abs(x - median_x))
    if mad == 0:
        return np.zeros_like(x)
    return 0.6745 * (x - median_x) / mad


def local_outliers(
    adata,
    metric: str = "detected",
    direction: str = "lower",
    n_neighbors: int = 36,
    samples: str = "sample_id",
    log: bool = True,
    cutoff: float = 3.0,
    knn_indices: Optional[np.ndarray] = None,
) -> "anndata.AnnData":

    """Detect local outliers using kNN-based modified z-scores.

    For each spot, computes a modified z-score of the metric within its k
    nearest neighbors using the MAD-based estimator (matching
    ``spatialEco::outliers``). Spots whose z-score exceeds the cutoff are
    flagged as outliers.

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data with spatial coordinates in
        ``adata.obsm['spatial']`` and QC metrics in ``adata.obs``.
    metric : str
        Column name in ``adata.obs`` for outlier detection.
    direction : str
        Direction of outlier detection: 'lower', 'higher', or 'both'.
    n_neighbors : int
        Number of nearest neighbors.
    samples : str
        Column name in ``adata.obs`` for sample IDs.
    log : bool
        Whether to log1p-transform the metric.
    cutoff : float
        Z-score cutoff for outlier detection (default 3).

    Returns
    -------
    anndata.AnnData
        Modified AnnData with ``{metric}_z`` (z-scores) and
        ``{metric}_outliers`` (boolean flags) added to ``adata.obs``.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Metric '{metric}' not found in adata.obs.")

    if direction not in ("lower", "higher", "both"):
        raise ValueError("'direction' must be 'lower', 'higher', or 'both'.")

    metric_to_use = metric
    if log:
        metric_log = f"{metric}_log"
        adata.obs[metric_log] = np.log1p(adata.obs[metric].values)
        metric_to_use = metric_log

    coords = adata.obsm["spatial"]
    values = adata.obs[metric_to_use].values.copy()
    sample_ids = adata.obs[samples].values

    z_scores = np.zeros(len(adata))
    outlier_flags = np.zeros(len(adata), dtype=bool)

    for sample in np.unique(sample_ids):
        mask = sample_ids == sample
        sample_coords = coords[mask]
        sample_values = values[mask]

        if knn_indices is not None:
            current_knn = knn_indices[mask] if len(knn_indices) == len(adata) else knn_indices
        else:
            current_knn = find_knn(sample_coords, n_neighbors)

        n_spots = len(sample_coords)
        sample_z = np.zeros(n_spots)

        for i in range(n_spots):
            neighbor_indices = current_knn[i]
            # R code uses neighbors only (NOT the focal spot)
            # spatialEco::outliers returns z-scores for all neighbor values
            # and [1] takes the first neighbor's z-score
            neighborhood = sample_values[neighbor_indices]
            z = _modified_zscore(neighborhood)
            sample_z[i] = z[0]  # z-score of the first neighbor

        # Handle non-finite values
        sample_z[~np.isfinite(sample_z)] = 0

        # Apply cutoff
        if direction == "higher":
            sample_outliers = sample_z > cutoff
        elif direction == "lower":
            sample_outliers = sample_z < -cutoff
        else:  # both
            sample_outliers = (sample_z > cutoff) | (sample_z < -cutoff)

        z_scores[mask] = sample_z
        outlier_flags[mask] = sample_outliers

    metric_z = f"{metric}_z"
    metric_outliers = f"{metric}_outliers"
    adata.obs[metric_z] = z_scores
    adata.obs[metric_outliers] = outlier_flags

    return adata
