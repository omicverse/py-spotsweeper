"""Shared utility functions for SpotSweeper."""

import numpy as np
from scipy.spatial import cKDTree


def find_knn(coords: np.ndarray, k: int) -> np.ndarray:
    """Find k nearest neighbors for each point.

    Parameters
    ----------
    coords : np.ndarray
        Spatial coordinates, shape (n_points, 2).
    k : int
        Number of nearest neighbors.

    Returns
    -------
    np.ndarray
        Index matrix of shape (n_points, k). Each row contains the indices
        of the k nearest neighbors (excluding the point itself).
    """
    tree = cKDTree(coords)
    # k+1 because the first neighbor is the point itself
    distances, indices = tree.query(coords, k=k + 1)
    # Remove the self-match (first column)
    return indices[:, 1:]
