# Benchmark R SpotSweeper functions
library(SpotSweeper)
library(SpatialExperiment)
library(STexampleData)
library(BiocNeighbors)
library(jsonlite)

cat("=== R Benchmark ===\n")

# Load data
spe <- STexampleData::Visium_humanDLPFC()
rownames(spe) <- rowData(spe)$gene_name
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))

# Export kNN for Python to use (to ensure same neighbors)
coords <- spatialCoords(spe)
dnn <- findKNN(coords, k = 36, warn.ties = FALSE)$index
write.csv(dnn, "data/r_knn_indices_k36.csv", row.names=FALSE)

results <- list()

# Benchmark localVariance
set.seed(42)
t1 <- system.time({
    spe_bench <- localVariance(spe, metric="subsets_Mito_percent", n_neighbors=36, name="lv", workers=1)
})
results$localVariance <- t1["elapsed"]
cat("localVariance:", t1["elapsed"], "s\n")

# Benchmark localOutliers
set.seed(42)
t2 <- system.time({
    spe_bench <- localOutliers(spe, metric="sum", direction="lower", log=TRUE, n_neighbors=36, workers=1)
})
results$localOutliers <- t2["elapsed"]
cat("localOutliers:", t2["elapsed"], "s\n")

# Benchmark flagVisiumOutliers
t3 <- system.time({
    spe_bench <- flagVisiumOutliers(spe)
})
results$flagVisiumOutliers <- t3["elapsed"]
cat("flagVisiumOutliers:", t3["elapsed"], "s\n")

# Benchmark findArtifacts (using DLPFC_artifact)
data(DLPFC_artifact, package="SpotSweeper")
spe_art <- DLPFC_artifact

# Export artifact kNN
art_coords <- spatialCoords(spe_art)
for (k in c(6, 18)) {
    dnn_art <- findKNN(art_coords, k = k, warn.ties = FALSE)$index
    write.csv(dnn_art, paste0("data/r_artifact_knn_k", k, ".csv"), row.names=FALSE)
}

set.seed(42)
t4 <- system.time({
    spe_art_bench <- findArtifacts(spe_art, mito_percent="expr_chrM_ratio", mito_sum="expr_chrM", n_order=2)
})
results$findArtifacts <- t4["elapsed"]
cat("findArtifacts:", t4["elapsed"], "s\n")

# Save results
write_json(results, "data/r_benchmark_times.json", auto_unbox=TRUE)
cat("\nResults saved to data/r_benchmark_times.json\n")

# Also save correlation reference outputs
set.seed(42)
spe <- localVariance(spe, metric="subsets_Mito_percent", n_neighbors=36, name="local_mito_variance_k36", workers=1)
spe <- localOutliers(spe, metric="sum", direction="lower", log=TRUE, n_neighbors=36, workers=1)
spe <- flagVisiumOutliers(spe)

ref <- list(
    local_variance_residuals = as.numeric(spe$local_mito_variance_k36),
    local_outlier_zscores = as.numeric(spe$sum_z),
    local_outlier_flags = as.logical(spe$sum_outliers),
    systematic_outlier_flags = as.logical(spe$systematic_outliers)
)
write_json(ref, "data/reference_outputs.json", auto_unbox=TRUE, digits=22)

set.seed(42)
spe_art <- findArtifacts(spe_art, mito_percent="expr_chrM_ratio", mito_sum="expr_chrM", n_order=2)
art_ref <- list(artifact_labels = as.logical(spe_art$artifact))
write_json(art_ref, "data/reference_artifact_outputs.json", auto_unbox=TRUE, digits=22)

cat("=== Done ===\n")
