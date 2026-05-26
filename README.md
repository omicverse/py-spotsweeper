# py-SpotSweeper

Spatially-aware quality control for spatial transcriptomics.

Python port of R/Bioconductor package [SpotSweeper](https://github.com/MicTott/SpotSweeper) v1.5.0.

## Install

```bash
pip install py-spotsweeper
```

## Quick start

```python
import anndata as ad
from spotsweeper import local_variance, local_outliers, find_artifacts, flag_visium_outliers

# Load your spatial data as AnnData with adata.obsm['spatial']
adata = ad.read_h5ad("your_data.h5ad")

# 1. Compute local variance of mitochondrial percentage
adata = local_variance(adata, metric="subsets_Mito_percent", n_neighbors=36)

# 2. Detect local outliers in library size
adata = local_outliers(adata, metric="sum", direction="lower", log=True)

# 3. Flag systematic Visium outliers
adata = flag_visium_outliers(adata)

# 4. Find artifacts (single sample only)
adata = find_artifacts(adata, mito_percent="subsets_Mito_percent", n_order=5)
```

## Functions

| Function | Description |
|---|---|
| `local_variance()` | Local variance of QC metrics using kNN + robust regression |
| `local_outliers()` | Outlier detection using MAD-based modified z-scores |
| `find_artifacts()` | Artifact detection via multi-scale variance + PCA + k-means |
| `flag_visium_outliers()` | Flag known systematic Visium outlier spots |
| `plot_qc_metrics()` | Spatial scatter plot of QC metrics |
| `plot_qc_pdf()` | Multi-page PDF of QC plots per sample |

## Speed vs R

| Function | R | Python | Speedup |
|---|---|---|---|
| localVariance | 8.9s | 0.5s | **19x** |
| localOutliers | 1.3s | 0.2s | **8x** |
| findArtifacts | 16s | 5.5s | **3x** |
| flagVisiumOutliers | 0.05s | 0.003s | **14x** |

## Parity

All functions pass parity gates against R reference outputs:

| Function | Metric | Value | Gate |
|---|---|---|---|
| localVariance | max abs err | 2.67e-6 | 1e-5 |
| localOutliers | F1 | 1.0 | 0.95 |
| findArtifacts | ARI | 1.0 | 0.95 |
| flagVisiumOutliers | exact match | 1.0 | 1.0 |

## Requirements

- Python >= 3.9
- numpy, scipy, pandas, anndata, scikit-learn, matplotlib

## License

MIT

## References

Totty et al. (2025) "SpotSweeper: spatially-aware quality control for spatial transcriptomics." Bioconductor.
