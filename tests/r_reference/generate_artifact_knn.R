# Generate kNN indices for DLPFC_artifact fixture
library(SpatialExperiment)
library(BiocNeighbors)

data(DLPFC_artifact, package="SpotSweeper")
spe <- DLPFC_artifact
coords <- spatialCoords(spe)

cat("Artifact data:", nrow(spe), "genes x", ncol(spe), "spots\n")

# Generate kNN for each order used in findArtifacts (n_order=2, hexagonal)
# Order 1: k = 3*1*2 = 6
# Order 2: k = 3*2*3 = 18

for (k in c(6, 18)) {
    dnn <- findKNN(coords, k = k, warn.ties = FALSE)$index
    fname <- paste0("data/r_artifact_knn_k", k, ".csv")
    write.csv(dnn, fname, row.names=FALSE)
    cat("Exported kNN for k=", k, "shape:", dim(dnn), "\n")
}

cat("Done\n")
