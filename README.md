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

### `local_variance()`

Calculate local variance of a QC metric within kNN neighborhoods, then regress out mean-variance bias using Huber robust linear regression (M-estimator). Returns adjusted residuals as local variance scores.

```python
from spotsweeper import local_variance

adata = local_variance(
    adata,
    metric="subsets_Mito_percent",  # column in adata.obs
    n_neighbors=36,                  # kNN size
    samples="sample_id",            # sample grouping column
    log=False,                       # log1p transform before computation
    name=None,                       # output column name, default: "{metric}_var"
    knn_indices=None,                # optional precomputed kNN index array
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | — | Spatial data with `adata.obsm['spatial']` and QC metrics in `adata.obs` |
| `metric` | str | `"expr_chrM_ratio"` | Column name in `adata.obs` to compute local variance for |
| `n_neighbors` | int | `36` | Number of nearest neighbors |
| `samples` | str | `"sample_id"` | Column name in `adata.obs` for sample IDs |
| `log` | bool | `False` | Whether to log1p-transform the metric before computation |
| `name` | str \| None | `None` | Output column name. Defaults to `"{metric}_var"` |
| `knn_indices` | np.ndarray \| None | `None` | Precomputed kNN indices (shape: n_spots x k). If None, computed internally |

**Output columns in `adata.obs`:** `{metric}_var` (or custom `name`)

---

### `local_outliers()`

Detect local outliers using MAD-based modified z-scores within kNN neighborhoods. For each spot, computes a modified z-score of the metric among its neighbors (matching R's `spatialEco::outliers`). Spots exceeding the z-score cutoff are flagged.

```python
from spotsweeper import local_outliers

adata = local_outliers(
    adata,
    metric="detected",       # column in adata.obs
    direction="lower",       # 'lower', 'higher', or 'both'
    n_neighbors=36,
    samples="sample_id",
    log=True,
    cutoff=3.0,
    knn_indices=None,
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | — | Spatial data with `adata.obsm['spatial']` |
| `metric` | str | `"detected"` | Column name in `adata.obs` for outlier detection |
| `direction` | str | `"lower"` | Direction: `"lower"`, `"higher"`, or `"both"` |
| `n_neighbors` | int | `36` | Number of nearest neighbors |
| `samples` | str | `"sample_id"` | Column name in `adata.obs` for sample IDs |
| `log` | bool | `True` | Whether to log1p-transform the metric |
| `cutoff` | float | `3.0` | Modified z-score cutoff for flagging outliers |
| `knn_indices` | np.ndarray \| None | `None` | Precomputed kNN indices |

**Output columns in `adata.obs`:**
- `{metric}_z` — modified z-scores
- `{metric}_outliers` — boolean outlier flags

---

### `find_artifacts()`

Identify regional artifacts by computing local mitochondrial variance at multiple neighborhood scales, running PCA on the variance matrix, and k-means clustering (k=2). The cluster with lower average local mito variance is labeled as artifact. **Input must contain exactly one sample.**

```python
from spotsweeper import find_artifacts

adata = find_artifacts(
    adata,
    mito_percent="subsets_Mito_percent",
    mito_sum="subsets_Mito_sum",
    samples="sample_id",
    n_order=5,
    shape="hexagonal",
    log=True,
    name="artifact",
    var_output=True,
    seed=42,
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | — | Spatial data (single sample only) |
| `mito_percent` | str | `"expr_chrM_ratio"` | Column name for mitochondrial percentage |
| `mito_sum` | str | `"expr_chrM"` | Column name for mitochondrial sum |
| `samples` | str | `"sample_id"` | Column name for sample IDs |
| `n_order` | int | `5` | Number of neighborhood scales. Each order `i` uses `3*i*(i+1)` neighbors (hex) or `4*i*(i+1)` (square) |
| `shape` | str | `"hexagonal"` | Neighborhood shape: `"hexagonal"` or `"square"` |
| `log` | bool | `True` | Whether to log1p-transform `mito_percent` |
| `name` | str | `"artifact"` | Output column name |
| `var_output` | bool | `True` | If True, keep intermediate variance columns (e.g. `k6`, `k12`, ...) |
| `seed` | int | `42` | Random seed for PCA and k-means |
| `knn_indices_dict` | dict \| None | `None` | Dict mapping `n_neighbors` to precomputed kNN arrays |

**Output columns in `adata.obs`:** `artifact` (boolean), optionally `k6`, `k12`, `k18`, ... (per-scale variances)

**Output in `adata.obsm`:** `X_pca_artifacts` (PCA embedding)

---

### `flag_visium_outliers()`

Flag known systematic Visium technical outlier spots by matching array coordinates against a built-in list of biased spots.

```python
from spotsweeper import flag_visium_outliers

adata = flag_visium_outliers(adata)
```

**Required columns in `adata.obs`:** `array_row`, `array_col`

**Output columns in `adata.obs`:** `systematic_outliers` (boolean)

---

### `plot_qc_metrics()`

Spatial scatter plot of a single QC metric, with optional outlier highlighting.

```python
from spotsweeper.plot_qc import plot_qc_metrics

fig = plot_qc_metrics(
    adata,
    metric="detected",
    outliers="detected_outliers",  # optional boolean column
    sample_id="sample_id",
    sample=None,                   # first sample by default
    point_size=2.0,
    colors=("white", "black"),
    stroke=1.0,
    ax=None,
)
```

---

### `plot_qc_pdf()`

Generate a multi-page PDF with one QC plot per sample.

```python
from spotsweeper.plot_qc import plot_qc_pdf

plot_qc_pdf(
    adata,
    metric="detected",
    outliers="detected_outliers",
    sample_id="sample_id",
    fname="qc_plots.pdf",
    point_size=2.0,
    colors=("white", "black"),
    stroke=1.0,
    width=5.0,
    height=5.0,
)
```

---

## R ↔ Python dictionary mapping

### Function mapping

| R (SpotSweeper) | Python (py-SpotSweeper) |
|---|---|
| `localVariance()` | `local_variance()` |
| `localOutliers()` | `local_outliers()` |
| `findArtifacts()` | `find_artifacts()` |
| `flagVisiumOutliers()` | `flag_visium_outliers()` |
| `plotQCmetrics()` | `plot_qc_metrics()` |
| `plotQCpdf()` | `plot_qc_pdf()` |

### Parameter mapping

#### `localVariance()` → `local_variance()`

| R parameter | Python parameter | Default (R / Py) | Notes |
|---|---|---|---|
| `spe` | `adata` | — | R: SpatialExperiment; Py: AnnData |
| `metric` | `metric` | `"expr_chrM_ratio"` / `"expr_chrM_ratio"` | — |
| `n_neighbors` | `n_neighbors` | `36` / `36` | — |
| `samples` | `samples` | `"sample_id"` / `"sample_id"` | — |
| `log` | `log` | `FALSE` / `False` | — |
| `name` | `name` | `NULL` / `None` | Default output: `"{metric}_var"` |
| `workers` | *(not supported)* | `1` / — | R uses BiocParallel; Py is single-threaded |
| — | `knn_indices` | — / `None` | Py-only: pass precomputed kNN to skip recomputation |

#### `localOutliers()` → `local_outliers()`

| R parameter | Python parameter | Default (R / Py) | Notes |
|---|---|---|---|
| `spe` | `adata` | — | R: SpatialExperiment/Seurat; Py: AnnData |
| `metric` | `metric` | `"detected"` / `"detected"` | — |
| `direction` | `direction` | `"lower"` / `"lower"` | `"lower"`, `"higher"`, `"both"` |
| `n_neighbors` | `n_neighbors` | `36` / `36` | — |
| `samples` | `samples` | `"sample_id"` / `"sample_id"` | — |
| `log` | `log` | `TRUE` / `True` | — |
| `cutoff` | `cutoff` | `3` / `3.0` | — |
| `workers` | *(not supported)* | `1` / — | R uses BiocParallel |
| `coords` | `knn_indices` | `NULL` / `None` | R: raw coord matrix; Py: precomputed kNN index array |

#### `findArtifacts()` → `find_artifacts()`

| R parameter | Python parameter | Default (R / Py) | Notes |
|---|---|---|---|
| `spe` | `adata` | — | Both require single-sample input |
| `mito_percent` | `mito_percent` | `"expr_chrM_ratio"` / `"expr_chrM_ratio"` | — |
| `mito_sum` | `mito_sum` | `"expr_chrM"` / `"expr_chrM"` | — |
| `samples` | `samples` | `"sample_id"` / `"sample_id"` | — |
| `n_order` | `n_order` | `5` / `5` | — |
| `shape` | `shape` | `"hexagonal"` / `"hexagonal"` | — |
| `log` | `log` | `TRUE` / `True` | — |
| `name` | `name` | `"artifact"` / `"artifact"` | — |
| `var_output` | `var_output` | `TRUE` / `True` | — |
| — | `seed` | — / `42` | Py-only: random seed for PCA/k-means |
| — | `knn_indices_dict` | — / `None` | Py-only: dict of precomputed kNN arrays |

#### `flagVisiumOutliers()` → `flag_visium_outliers()`

| R parameter | Python parameter | Notes |
|---|---|---|
| `spe` | `adata` | R: SpatialExperiment; Py: AnnData |

No additional parameters in either version. Requires `array_row` and `array_col` columns.

#### `plotQCmetrics()` → `plot_qc_metrics()`

| R parameter | Python parameter | Default (R / Py) | Notes |
|---|---|---|---|
| `spe` | `adata` | — | — |
| `sample_id` | `sample_id` | `"sample_id"` / `"sample_id"` | — |
| `sample` | `sample` | first sample / `None` | — |
| `metric` | `metric` | `"detected"` / `"detected"` | — |
| `outliers` | `outliers` | `NULL` / `None` | — |
| `point_size` | `point_size` | `2` / `2.0` | — |
| `colors` | `colors` | `c("white","black")` / `("white","black")` | — |
| `stroke` | `stroke` | `1` / `1.0` | — |
| — | `ax` | — / `None` | Py-only: pass matplotlib Axes |

#### `plotQCpdf()` → `plot_qc_pdf()`

| R parameter | Python parameter | Default (R / Py) | Notes |
|---|---|---|---|
| `spe` | `adata` | — | — |
| `sample_id` | `sample_id` | `"sample_id"` / `"sample_id"` | — |
| `metric` | `metric` | `"detected"` / `"detected"` | — |
| `outliers` | `outliers` | `"local_outliers"` / `None` | — |
| `colors` | `colors` | `c("white","black")` / `("white","black")` | — |
| `stroke` | `stroke` | `1` / `1.0` | — |
| `point_size` | `point_size` | `2` / `2.0` | — |
| `width` | `width` | `5` / `5.0` | — |
| `height` | `height` | `5` / `5.0` | — |
| `fname` | `fname` | — / `"qc_plots.pdf"` | — |

### Object model mapping

| R | Python |
|---|---|
| `SpatialExperiment` | `AnnData` |
| `colData(spe)` | `adata.obs` |
| `spatialCoords(spe)` | `adata.obsm['spatial']` |
| `reducedDim(spe, "PCA")` | `adata.obsm['X_pca_artifacts']` |
| `BiocNeighbors::findKNN()` | `scipy.spatial.cKDTree` |
| `MASS::rlm()` (Huber) | `_huber_rlm()` (IRLS, matching R) |
| `spatialEco::outliers()` | `_modified_zscore()` (MAD-based) |

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
- numpy, scipy, pandas, anndata, scikit-learn, statsmodels, matplotlib

## License

MIT

## References

Totty et al. (2025) "SpotSweeper: spatially-aware quality control for spatial transcriptomics." Bioconductor.
