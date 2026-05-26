# Generate R reference data with HIGH PRECISION + kNN indices
library(SpotSweeper)
library(SpatialExperiment)
library(STexampleData)
library(BiocNeighbors)
library(jsonlite)

cat("=== Generating high-precision reference data ===\n")

# Load and prepare data
spe <- STexampleData::Visium_humanDLPFC()
rownames(spe) <- rowData(spe)$gene_name
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))

cat("Data:", nrow(spe), "genes x", ncol(spe), "spots\n")

# Export kNN indices for verification
coords <- spatialCoords(spe)
dnn <- findKNN(coords, k = 36, warn.ties = FALSE)$index
cat("kNN indices shape:", dim(dnn), "\n")
cat("kNN[1,1:5]:", dnn[1, 1:5], "\n")

# Export kNN as CSV (1-based indices)
write.csv(dnn, "data/r_knn_indices_k36.csv", row.names=FALSE)

# Run localVariance
set.seed(42)
spe <- localVariance(spe,
    metric = "subsets_Mito_percent",
    n_neighbors = 36,
    name = "local_mito_variance_k36",
    workers = 1
)

# Run localOutliers
spe <- localOutliers(spe,
    metric = "sum",
    direction = "lower",
    log = TRUE,
    n_neighbors = 36,
    workers = 1
)

# Run flagVisiumOutliers
spe <- flagVisiumOutliers(spe)

# Export with FULL precision
reference <- list(
    local_variance_residuals = as.numeric(spe$local_mito_variance_k36),
    local_outlier_zscores = as.numeric(spe$sum_z),
    local_outlier_flags = as.logical(spe$sum_outliers),
    systematic_outlier_flags = as.logical(spe$systematic_outliers),
    local_outlier_log_values = as.numeric(spe$sum_log)
)

json_text <- toJSON(reference, auto_unbox=TRUE, digits=22, pretty=FALSE)
writeLines(json_text, "data/reference_outputs.json")
cat("High-precision reference exported.\n")

# Artifact reference
data(DLPFC_artifact, package="SpotSweeper")
spe_artifact <- DLPFC_artifact
set.seed(42)
spe_artifact <- findArtifacts(spe_artifact,
    mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM",
    n_order = 2,
    name = "artifact"
)

artifact_ref <- list(
    artifact_labels = as.logical(spe_artifact$artifact)
)
json_text2 <- toJSON(artifact_ref, auto_unbox=TRUE, digits=22, pretty=FALSE)
writeLines(json_text2, "data/reference_artifact_outputs.json")
cat("Artifact reference exported.\n")
cat("=== Done ===\n")
