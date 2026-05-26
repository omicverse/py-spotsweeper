"""Local variance computation using kNN neighborhoods."""

from typing import Optional
import numpy as np
from spotsweeper.utils import find_knn


def _huber_rlm(y, X, c=1.345, max_iter=20, tol=1e-5):
    """Robust linear regression matching R's MASS::rlm with Huber psi.

    Implements IRLS (Iteratively Reweighted Least Squares) with the Huber
    weight function and MAD scale estimation, exactly matching MASS::rlm.

    Parameters
    ----------
    y : np.ndarray
        Response variable.
    X : np.ndarray
        Design matrix (should include intercept column).
    c : float
        Tuning constant for Huber function (default 1.345, matching R).
    max_iter : int
        Maximum IRLS iterations (default 20, matching R).
    tol : float
        Convergence tolerance (default 1e-5, matching R).

    Returns
    -------
    np.ndarray
        Raw residuals (y - X @ beta).
    """
    # Initial OLS fit
    try:
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        beta = np.zeros(X.shape[1])

    for iteration in range(max_iter):
        resid = y - X @ beta

        # MAD scale estimate (matching R's MASS::rlm)
        s = 1.4826 * np.median(np.abs(resid))
        if s < 1e-10:
            break

        # Huber weights: w = psi(r/s) / (r/s)
        u = resid / s
        abs_u = np.abs(u)
        weights = np.where(abs_u <= c, 1.0, c / abs_u)

        # Weighted least squares: solve X^T W X beta = X^T W y
        try:
            W = np.diag(weights)
            XtWX = X.T @ W @ X
            XtWy = X.T @ W @ y
            beta_new = np.linalg.lstsq(XtWX, XtWy, rcond=None)[0]
        except np.linalg.LinAlgError:
            break

        # Check convergence
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break

        beta = beta_new

    return y - X @ beta


def local_variance(
    adata,
    metric: str = "expr_chrM_ratio",
    n_neighbors: int = 36,
    samples: str = "sample_id",
    log: bool = False,
    name: Optional[str] = None,
    knn_indices: Optional[np.ndarray] = None,
) -> "anndata.AnnData":
    """Calculate local variance based on kNN neighborhoods.

    For each spot, computes the variance of a metric within its k nearest
    neighbors, then regresses out mean-variance bias using robust linear
    regression (Huber M-estimator). The residuals represent local variance
    adjusted for the mean-variance relationship.

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data with spatial coordinates in
        ``adata.obsm['spatial']`` and QC metrics in ``adata.obs``.
    metric : str
        Column name in ``adata.obs`` to compute local variance for.
    n_neighbors : int
        Number of nearest neighbors.
    samples : str
        Column name in ``adata.obs`` for sample IDs.
    log : bool
        Whether to log1p-transform the metric before computation.
    name : str or None
        Name for the output column. Defaults to ``"{metric}_var"``.

    Returns
    -------
    anndata.AnnData
        Modified AnnData with local variance column added to ``adata.obs``.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Metric '{metric}' not found in adata.obs.")

    if samples not in adata.obs.columns:
        raise ValueError(f"Samples column '{samples}' not found in adata.obs.")

    if name is None:
        name = f"{metric}_var"

    metric_to_use = metric
    if log:
        metric_log = f"{metric}_log"
        adata.obs[metric_log] = np.log1p(adata.obs[metric].values)
        metric_to_use = metric_log

    coords = adata.obsm["spatial"]
    values = adata.obs[metric_to_use].values.copy()
    sample_ids = adata.obs[samples].values
    residuals = np.full(len(adata), np.nan)

    for sample_idx, sample in enumerate(np.unique(sample_ids)):
        mask = sample_ids == sample
        sample_coords = coords[mask]
        sample_values = values[mask]

        # Find kNN (use provided indices or compute)
        if knn_indices is not None:
            if knn_indices.ndim == 1:
                knn_indices = knn_indices.reshape(1, -1)
            current_knn = knn_indices[mask] if len(knn_indices) == len(adata) else knn_indices
        else:
            current_knn = find_knn(sample_coords, n_neighbors)

        # Compute local variance and mean for each spot
        n_spots = len(sample_coords)
        var_vals = np.zeros(n_spots)
        mean_vals = np.zeros(n_spots)

        for i in range(n_spots):
            neighbor_indices = current_knn[i]
            # Include the focal spot itself (matching R behavior)
            all_indices = np.concatenate([[i], neighbor_indices])
            neighborhood = sample_values[all_indices]
            var_vals[i] = np.var(neighborhood, ddof=1)
            mean_vals[i] = np.mean(neighborhood)

        # Handle non-finite values (matching R: stats_matrix[!is.finite(stats_matrix)] <- 0)
        var_vals[~np.isfinite(var_vals)] = 0
        mean_vals[~np.isfinite(mean_vals)] = 0

        # Robust linear regression: log2(var) ~ mean
        # R's formula: rlm(log2(var) ~ mean)
        # log2(0) = -Inf in both R and Python
        log2_var = np.log2(var_vals)

        X = np.column_stack([np.ones(n_spots), mean_vals])
        resid_vals = _huber_rlm(log2_var, X)

        # Direct assignment to the masked slice of the original array
        residuals[np.where(mask)[0]] = resid_vals

    adata.obs[name] = residuals

    return adata
