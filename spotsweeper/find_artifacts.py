"""Artifact detection using local mito variance and kmeans clustering."""

from typing import Dict, Optional
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from spotsweeper.local_variance import local_variance


def find_artifacts(
    adata,
    mito_percent: str = "expr_chrM_ratio",
    mito_sum: str = "expr_chrM",
    samples: str = "sample_id",
    n_order: int = 5,
    shape: str = "hexagonal",
    log: bool = True,
    name: str = "artifact",
    var_output: bool = True,
    seed: int = 42,
    knn_indices_dict: Optional[Dict[int, np.ndarray]] = None,
) -> "anndata.AnnData":
    """Identify and annotate artifacts in spatial transcriptomics data.

    Computes local mitochondrial variance at multiple neighborhood scales,
    performs PCA on the variance matrix, and uses k-means clustering to
    identify artifact spots (those with consistently low local variance).

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data. Must contain exactly one sample.
    mito_percent : str
        Column name for mitochondrial percentage.
    mito_sum : str
        Column name for mitochondrial sum.
    samples : str
        Column name for sample IDs.
    n_order : int
        Number of neighborhood orders for local variance calculation.
    shape : str
        Neighborhood shape: 'hexagonal' or 'square'.
    log : bool
        Whether to log1p-transform the mito_percent.
    name : str
        Name for the artifact column in ``adata.obs``.
    var_output : bool
        Whether to include local variances in the output.
    seed : int
        Random seed for kmeans clustering.

    Returns
    -------
    anndata.AnnData
        Modified AnnData with artifact annotations in ``adata.obs``.
    """
    if mito_percent not in adata.obs.columns:
        raise ValueError(f"'{mito_percent}' not found in adata.obs.")
    if mito_sum not in adata.obs.columns:
        raise ValueError(f"'{mito_sum}' not found in adata.obs.")
    if samples not in adata.obs.columns:
        raise ValueError(f"'{samples}' not found in adata.obs.")
    if shape not in ("hexagonal", "square"):
        raise ValueError("'shape' must be 'hexagonal' or 'square'.")

    unique_samples = adata.obs[samples].unique()
    if len(unique_samples) > 1:
        raise ValueError("Input must contain only one sample for find_artifacts.")

    # Compute local mito variance at each order
    var_matrix = np.zeros((len(adata), n_order))

    for i in range(1, n_order + 1):
        if shape == "hexagonal":
            n_neighbors = 3 * i * (i + 1)
        else:
            n_neighbors = 4 * i * (i + 1)

        tmp_name = f"k{n_neighbors}"

        # Compute local variance for this order
        knn = knn_indices_dict.get(n_neighbors) if knn_indices_dict else None
        adata = local_variance(
            adata,
            metric=mito_percent,
            n_neighbors=n_neighbors,
            samples=samples,
            log=log,
            name=tmp_name,
            knn_indices=knn,
        )

        var_matrix[:, i - 1] = adata.obs[tmp_name].values

    # Build PCA input: var_matrix + mito_percent + mito_sum
    pca_input = np.column_stack([
        var_matrix,
        adata.obs[mito_percent].values,
        adata.obs[mito_sum].values,
    ])

    # Run PCA with scaling (matching R's prcomp(center=TRUE, scale.=TRUE))
    pca = PCA(n_components=pca_input.shape[1], random_state=seed)
    # Standardize each column to mean=0, std=1 (matching R's scale.)
    pca_mean = pca_input.mean(axis=0)
    pca_std = pca_input.std(axis=0, ddof=1)
    pca_std[pca_std == 0] = 1  # avoid division by zero
    pca_input_scaled = (pca_input - pca_mean) / pca_std
    pca_result = pca.fit_transform(pca_input_scaled)
    adata.obsm["X_pca_artifacts"] = pca_result

    # K-means clustering (k=2)
    kmeans = KMeans(n_clusters=2, n_init=25, random_state=seed)
    cluster_labels = kmeans.fit_predict(pca_result)

    # Determine which cluster is the artifact (lower mean variance)
    # Use the k18 (hex) or k24 (square) variance column
    if shape == "hexagonal":
        ref_var_name = "k18"
    else:
        ref_var_name = "k24"

    ref_var = adata.obs[ref_var_name].values
    clus0_mean = np.mean(ref_var[cluster_labels == 0])
    clus1_mean = np.mean(ref_var[cluster_labels == 1])

    # Artifact cluster has lower average variance
    artifact_cluster = 0 if clus0_mean < clus1_mean else 1

    # Create artifact annotation
    adata.obs[name] = cluster_labels == artifact_cluster

    # Clean up intermediate columns if not requested
    if not var_output:
        for i in range(1, n_order + 1):
            if shape == "hexagonal":
                n_neighbors = 3 * i * (i + 1)
            else:
                n_neighbors = 4 * i * (i + 1)
            col_name = f"k{n_neighbors}"
            if col_name in adata.obs.columns:
                adata.obs.drop(columns=[col_name], inplace=True)

    return adata
