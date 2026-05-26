"""QC visualization functions for spatial transcriptomics data."""

from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def plot_qc_metrics(
    adata,
    metric: str = "detected",
    outliers: Optional[str] = None,
    sample_id: str = "sample_id",
    sample: Optional[str] = None,
    point_size: float = 2.0,
    colors: tuple[str, ...] = ("white", "black"),
    stroke: float = 1.0,
    ax: Optional[plt.Axes] = None,
) -> plt.Figure:
    """Plot QC metrics on spatial coordinates.

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data with spatial coordinates in
        ``adata.obsm['spatial']``.
    metric : str
        Column name in ``adata.obs`` to visualize.
    outliers : str or None
        Column name for outlier flags (boolean). If provided, outlier spots
        are highlighted with red borders.
    sample_id : str
        Column name for sample IDs.
    sample : str or None
        Specific sample to plot. If None, uses the first sample.
    point_size : float
        Size of scatter points.
    colors : tuple
        Color gradient endpoints (2 colors) or low/mid/high (3 colors).
    stroke : float
        Border width for outlier points.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates a new figure.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the plot.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Metric '{metric}' not found in adata.obs.")

    # Subset to sample
    if sample is None:
        sample = adata.obs[sample_id].unique()[0]

    mask = adata.obs[sample_id] == sample
    adata_sub = adata[mask]

    coords = adata_sub.obsm["spatial"]
    values = adata_sub.obs[metric].values

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    else:
        fig = ax.get_figure()

    # Create color map
    if len(colors) == 2:
        cmap = LinearSegmentedColormap.from_list("qc", colors)
    else:
        cmap = LinearSegmentedColormap.from_list("qc", colors)

    # Flip y-axis to match spatial orientation
    scatter = ax.scatter(
        coords[:, 0],
        -coords[:, 1],
        c=values,
        cmap=cmap,
        s=point_size,
        edgecolors="none",
    )

    # Add outlier highlights
    if outliers is not None and outliers in adata_sub.obs.columns:
        outlier_mask = adata_sub.obs[outliers].values.astype(bool)
        if outlier_mask.any():
            ax.scatter(
                coords[outlier_mask, 0],
                -coords[outlier_mask, 1],
                s=point_size * 1.5,
                facecolors="none",
                edgecolors="red",
                linewidths=stroke,
            )

    ax.set_title(f"Sample: {sample}")
    ax.set_aspect("equal")
    plt.colorbar(scatter, ax=ax, label=metric)

    return fig


def plot_qc_pdf(
    adata,
    metric: str = "detected",
    outliers: Optional[str] = None,
    sample_id: str = "sample_id",
    fname: str = "qc_plots.pdf",
    point_size: float = 2.0,
    colors: tuple[str, ...] = ("white", "black"),
    stroke: float = 1.0,
    width: float = 5.0,
    height: float = 5.0,
) -> None:
    """Generate a PDF with QC plots for each sample.

    Parameters
    ----------
    adata : anndata.AnnData
        Spatial transcriptomics data.
    metric : str
        Column name for the metric to visualize.
    outliers : str or None
        Column name for outlier flags.
    sample_id : str
        Column name for sample IDs.
    fname : str
        Output PDF file path.
    point_size : float
        Size of scatter points.
    colors : tuple
        Color gradient.
    stroke : float
        Border width for outlier points.
    width : float
        Figure width in inches.
    height : float
        Figure height in inches.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    unique_samples = adata.obs[sample_id].unique()

    with PdfPages(fname) as pdf:
        for sample in unique_samples:
            fig = plot_qc_metrics(
                adata,
                metric=metric,
                outliers=outliers,
                sample_id=sample_id,
                sample=sample,
                point_size=point_size,
                colors=colors,
                stroke=stroke,
            )
            fig.set_size_inches(width, height)
            pdf.savefig(fig)
            plt.close(fig)
