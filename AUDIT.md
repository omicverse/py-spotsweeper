# AUDIT.md — R Function Coverage

## Exported R Functions

| # | R Function | Ported | Python Function | Notes |
|---|-----------|--------|-----------------|-------|
| 1 | `findArtifacts` | Yes | `find_artifacts()` | Full port. Python adds `seed` and `knn_indices_dict` params. |
| 2 | `flagVisiumOutliers` | Yes | `flag_visium_outliers()` | Full port. Uses CSV instead of R .rds. |
| 3 | `getMetadata` | Skipped | — | R-only Seurat/SpatialExperiment shim. AnnData uses `adata.obs`. |
| 4 | `getSpatialCoords` | Skipped | — | R-only shim. AnnData uses `adata.obsm['spatial']`. |
| 5 | `localOutliers` | Yes | `local_outliers()` | Full port. Reimplements `spatialEco::outliers`. |
| 6 | `localVariance` | Yes | `local_variance()` | Full port. Reimplements `MASS::rlm`. |
| 7 | `plotQCmetrics` | Yes | `plot_qc_metrics()` | Full port. matplotlib instead of escheR/ggplot2. |
| 8 | `plotQCpdf` | Yes | `plot_qc_pdf()` | Full port. PdfPages instead of R pdf()/dev.off(). |
| 9 | `setMetadata` | Skipped | — | R-only shim. Direct `adata.obs` assignment. |
| 10 | `subsetSpatialObject` | Skipped | — | R-only shim. Direct AnnData slicing. |
| 11 | `validateMetadataColumns` | Skipped | — | R-only shim. Inline `raise ValueError()`. |

## Summary

- **6/11 ported** (55%) — all algorithmic functions
- **5/11 skipped** (45%) — all Seurat/SpatialExperiment compatibility shims
- **0 unported** core functions

## Rationale for skipped functions

All 5 skipped functions (`getMetadata`, `getSpatialCoords`, `setMetadata`, `subsetSpatialObject`, `validateMetadataColumns`) exist in `R/seurat_compatibility.R` to provide a dual Seurat/SpatialExperiment abstraction. Python targets only AnnData, where these operations are native (`adata.obs`, `adata.obsm['spatial']`, `adata[mask]`). No abstraction layer needed.

## Audit classification

**Class A** — Pure translation. No algorithmic modifications, no acceleration rewrites.
