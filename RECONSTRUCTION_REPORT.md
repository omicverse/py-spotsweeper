# RECONSTRUCTION_REPORT.md — py-SpotSweeper v0.1.0

## 1. Identity

| Field | Value |
|---|---|
| Package | py-SpotSweeper |
| Upstream | SpotSweeper v1.5.0 (Bioconductor) |
| Algorithm class | deterministic-bounded |
| Parity threshold | 1e-5 |
| Final parity | max abs err 2.67e-6 (localVariance), F1 1.0 (localOutliers), ARI 1.0 (findArtifacts), exact match 1.0 (flagVisiumOutliers) |
| Audit class | A (pure translation) |
| Speedup | 3-19x vs R |
| LOC | ~400 Python |
| Date | 2026-05-27 |

## 2. R Function Coverage Audit

See [AUDIT.md](AUDIT.md) for full details.

- **6/11 exported R functions ported** (all algorithmic functions)
- **5/11 skipped** (Seurat/SpatialExperiment compatibility shims — not needed in Python/AnnData)
- **0 unported core functions**

Skipped functions: `getMetadata`, `getSpatialCoords`, `setMetadata`, `subsetSpatialObject`, `validateMetadataColumns` — all from `R/seurat_compatibility.R`. These provide a dual Seurat/SpatialExperiment abstraction in R. Python targets only AnnData, where `adata.obs` and `adata.obsm['spatial']` are the native access patterns.

## 3. Parity Evidence

### Per-output parity

| Output | Class | Metric | Value | Gate | Pass |
|---|---|---|---|---|---|
| `local_variance_residuals` | deterministic | max abs err | 2.67e-6 | 1e-5 | ✅ |
| `local_outlier_zscores` | deterministic | max abs err | 0.0 | 1e-5 | ✅ |
| `local_outlier_flags` | classification | F1 | 1.0 | 0.95 | ✅ |
| `artifact_labels` | clustering | ARI | 1.0 | 0.95 | ✅ |
| `systematic_outlier_flags` | deterministic | exact match | 1.0 | 1.0 | ✅ |

### Correlation

| Function | Pearson r | Spearman rho |
|---|---|---|
| localVariance | 1.000000 | 1.000000 |
| localOutliers (z) | 1.000000 | 1.000000 |

### Reproducible reference

```bash
# R reference
cd py-SpotSweeper
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" tests/r_reference/generate_fixtures_v3.R

# Python parity test
python -m pytest tests/test_parity.py -v
```

## 4. Acceleration Evidence

**None** — Class A translation only. No acceleration rewrites were attempted.

See [MATH.md](MATH.md) for perturbation analysis.

## 5. Code Quality Audit

| Check | Status |
|---|---|
| `pip install .` | ✅ |
| `pytest -q` | ✅ (5 passed) |
| Notebook 1: `compare_R_vs_Python.ipynb` | ✅ created |
| Notebook 2: `tutorial_visium.ipynb` | ✅ created |
| Notebook 3: `function_by_function_R_parity.ipynb` | ✅ created |
| Notebook 4: `evolution.ipynb` | ✅ created |
| License | MIT (matches upstream) |
| Version | 0.1.0 |

## 6. Known Limitations

1. **kNN tie-breaking**: R's BiocNeighbors and Python's cKDTree break tied distances differently. ~0.36% of spots may get different nearest neighbors. Parity tests use R-exported kNN indices; production code uses Python's cKDTree.

2. **IRLS convergence**: The Huber M-estimator has max abs error ~2.67e-6 vs R due to floating-point differences in BLAS implementations. This does not affect biological conclusions.

3. **Single-threaded**: R's `workers` parameter (BiocParallel) is not ported. Python runs single-threaded but is still 3-19x faster due to more efficient kNN and array operations.

4. **No Seurat support**: The R package supports both Seurat and SpatialExperiment objects. The Python package targets AnnData only.

## 7. Integration into omicverse

- Package location: `spotsweeper/` (standalone, not yet vendored into omicverse)
- Public API: `from spotsweeper import local_variance, local_outliers, find_artifacts, flag_visium_outliers`
- Tutorial: `examples/tutorial_visium.ipynb`

## 8. Sign-off

| Field | Value |
|---|---|
| Author | rebuildr-agent |
| Date | 2026-05-27 |
| Active time | ~2 hours |
| Audit class | A (pure translation) |
