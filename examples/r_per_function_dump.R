# R driver for Notebook 3 — per-function output dump
# Runs each SpotSweeper function individually and dumps outputs as JSON

library(SpotSweeper)
library(SpatialExperiment)
library(STexampleData)
library(BiocNeighbors)
library(jsonlite)

cat("=== Per-function R output dump ===\n")

# Load data
spe <- STexampleData::Visium_humanDLPFC()
rownames(spe) <- rowData(spe)$gene_name
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))

# Export kNN
coords <- spatialCoords(spe)
dnn <- findKNN(coords, k = 36, warn.ties = FALSE)$index
write.csv(dnn, "data/r_knn_indices_k36.csv", row.names=FALSE)

results <- list()

# 1. localVariance
cat("1. localVariance...\n")
set.seed(42)
spe_lv <- localVariance(spe, metric="subsets_Mito_percent", n_neighbors=36,
                        name="local_mito_variance_k36", workers=1)
results$localVariance <- list(
    residuals = as.numeric(colData(spe_lv)$local_mito_variance_k36)
)

# 2. localOutliers
cat("2. localOutliers...\n")
spe_lo <- localOutliers(spe, metric="sum", direction="lower", log=TRUE,
                        n_neighbors=36, workers=1)
results$localOutliers <- list(
    z_scores = as.numeric(colData(spe_lo)$sum_z),
    outlier_flags = as.logical(colData(spe_lo)$sum_outliers),
    log_values = as.numeric(colData(spe_lo)$sum_log)
)

# 3. findArtifacts
cat("3. findArtifacts...\n")
data(DLPFC_artifact, package="SpotSweeper")
spe_art <- DLPFC_artifact

# Export artifact kNN
art_coords <- spatialCoords(spe_art)
for (k in c(6, 18)) {
    dnn_art <- findKNN(art_coords, k = k, warn.ties = FALSE)$index
    write.csv(dnn_art, paste0("data/r_artifact_knn_k", k, ".csv"), row.names=FALSE)
}

set.seed(42)
spe_art <- findArtifacts(spe_art, mito_percent="expr_chrM_ratio",
                         mito_sum="expr_chrM", n_order=2)
results$findArtifacts <- list(
    artifact_labels = as.logical(colData(spe_art)$artifact)
)

# 4. flagVisiumOutliers
cat("4. flagVisiumOutliers...\n")
spe_fv <- flagVisiumOutliers(spe)
results$flagVisiumOutliers <- list(
    systematic_outliers = as.logical(colData(spe_fv)$systematic_outliers)
)

# Save
json_text <- toJSON(results, auto_unbox=TRUE, digits=22, pretty=TRUE)
writeLines(json_text, "data/r_per_function_outputs.json")
cat("Saved to data/r_per_function_outputs.json\n")
cat("=== Done ===\n")
