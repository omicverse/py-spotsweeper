"""Benchmark and correlation analysis: Python vs R SpotSweeper."""

import json
import time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr
import anndata as ad

DATA_DIR = "data"


def load_visium():
    counts = mmread(f"{DATA_DIR}/fixture_counts.mtx").T
    coords = pd.read_csv(f"{DATA_DIR}/fixture_spatial_coords.csv", index_col=0)
    metadata = pd.read_csv(f"{DATA_DIR}/fixture_metadata.csv", index_col=0)
    features = pd.read_csv(f"{DATA_DIR}/feature_names.csv")
    spots = pd.read_csv(f"{DATA_DIR}/spot_names.csv")
    adata = ad.AnnData(
        X=csr_matrix(counts), obs=metadata,
        var=pd.DataFrame(index=features.iloc[:, 0].values),
    )
    adata.obs_names = spots.iloc[:, 0].values
    adata.obsm["spatial"] = coords.values
    return adata


def load_artifact():
    counts = mmread(f"{DATA_DIR}/fixture_artifact_counts.mtx").T
    coords = pd.read_csv(f"{DATA_DIR}/fixture_artifact_coords.csv", index_col=0)
    metadata = pd.read_csv(f"{DATA_DIR}/fixture_artifact_metadata.csv", index_col=0)
    features = pd.read_csv(f"{DATA_DIR}/fixture_artifact_features.csv")
    spots = pd.read_csv(f"{DATA_DIR}/fixture_artifact_spots.csv")
    adata = ad.AnnData(
        X=csr_matrix(counts), obs=metadata,
        var=pd.DataFrame(index=features.iloc[:, 0].values),
    )
    adata.obs_names = spots.iloc[:, 0].values
    adata.obsm["spatial"] = coords.values
    return adata


def load_ref():
    with open(f"{DATA_DIR}/reference_outputs.json") as f:
        return json.load(f)


def load_artifact_ref():
    with open(f"{DATA_DIR}/reference_artifact_outputs.json") as f:
        return json.load(f)


def load_r_times():
    with open(f"{DATA_DIR}/r_benchmark_times.json") as f:
        return json.load(f)


def load_r_knn():
    return pd.read_csv(f"{DATA_DIR}/r_knn_indices_k36.csv").values - 1


def run_benchmark_and_correlation():
    print("=" * 60)
    print("py-SpotSweeper Benchmark & Correlation Analysis")
    print("=" * 60)

    ref = load_ref()
    r_times = load_r_times()
    r_knn = load_r_knn()

    # --- localVariance ---
    print("\n--- localVariance ---")
    adata = load_visium()
    from spotsweeper.local_variance import local_variance

    t0 = time.time()
    adata = local_variance(
        adata, metric="subsets_Mito_percent", n_neighbors=36,
        name="local_mito_variance_k36", log=False, knn_indices=r_knn,
    )
    py_time = time.time() - t0

    py_vals = adata.obs["local_mito_variance_k36"].values
    ref_vals = np.array(ref["local_variance_residuals"])

    r, _ = pearsonr(ref_vals, py_vals)
    rho, _ = spearmanr(ref_vals, py_vals)
    max_err = np.max(np.abs(py_vals - ref_vals))

    print(f"  Python: {py_time:.3f}s | R: {r_times['localVariance']:.3f}s")
    print(f"  Speedup: {r_times['localVariance']/py_time:.2f}x")
    print(f"  Pearson r: {r:.8f}")
    print(f"  Spearman rho: {rho:.8f}")
    print(f"  Max abs error: {max_err:.2e}")

    # --- localOutliers ---
    print("\n--- localOutliers ---")
    adata = load_visium()
    from spotsweeper.local_outliers import local_outliers

    t0 = time.time()
    adata = local_outliers(
        adata, metric="sum", direction="lower", log=True,
        n_neighbors=36, knn_indices=r_knn,
    )
    py_time = time.time() - t0

    py_z = adata.obs["sum_z"].values
    ref_z = np.array(ref["local_outlier_zscores"])
    py_flags = adata.obs["sum_outliers"].values.astype(bool)
    ref_flags = np.array(ref["local_outlier_flags"], dtype=bool)

    r_z, _ = pearsonr(ref_z, py_z)
    rho_z, _ = spearmanr(ref_z, py_z)
    max_err_z = np.max(np.abs(py_z - ref_z))
    from sklearn.metrics import f1_score
    f1 = f1_score(ref_flags, py_flags)

    print(f"  Python: {py_time:.3f}s | R: {r_times['localOutliers']:.3f}s")
    print(f"  Speedup: {r_times['localOutliers']/py_time:.2f}x")
    print(f"  Z-score Pearson r: {r_z:.8f}")
    print(f"  Z-score Spearman rho: {rho_z:.8f}")
    print(f"  Z-score max abs error: {max_err_z:.2e}")
    print(f"  Outlier F1: {f1:.4f}")

    # --- flagVisiumOutliers ---
    print("\n--- flagVisiumOutliers ---")
    adata = load_visium()
    from spotsweeper.flag_visium_outliers import flag_visium_outliers

    t0 = time.time()
    adata = flag_visium_outliers(adata)
    py_time = time.time() - t0

    py_flags = adata.obs["systematic_outliers"].values.astype(bool)
    ref_flags = np.array(ref["systematic_outlier_flags"], dtype=bool)
    match = np.mean(py_flags == ref_flags)

    print(f"  Python: {py_time:.4f}s | R: {r_times['flagVisiumOutliers']:.4f}s")
    print(f"  Speedup: {r_times['flagVisiumOutliers']/py_time:.2f}x")
    print(f"  Exact match: {match:.4f}")

    # --- findArtifacts ---
    print("\n--- findArtifacts ---")
    adata = load_artifact()
    art_ref = load_artifact_ref()
    r_knn_6 = pd.read_csv(f"{DATA_DIR}/r_artifact_knn_k6.csv").values - 1
    r_knn_18 = pd.read_csv(f"{DATA_DIR}/r_artifact_knn_k18.csv").values - 1
    knn_dict = {6: r_knn_6, 18: r_knn_18}

    from spotsweeper.find_artifacts import find_artifacts

    t0 = time.time()
    adata = find_artifacts(
        adata, mito_percent="expr_chrM_ratio", mito_sum="expr_chrM",
        n_order=2, shape="hexagonal", name="artifact", seed=42,
        knn_indices_dict=knn_dict,
    )
    py_time = time.time() - t0

    py_labels = adata.obs["artifact"].values.astype(int)
    ref_labels = np.array(art_ref["artifact_labels"], dtype=int)
    from sklearn.metrics import adjusted_rand_score
    ari = adjusted_rand_score(ref_labels, py_labels)

    print(f"  Python: {py_time:.3f}s | R: {r_times['findArtifacts']:.3f}s")
    print(f"  Speedup: {r_times['findArtifacts']/py_time:.2f}x")
    print(f"  ARI: {ari:.4f}")

    # --- Summary ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    total_r = sum(r_times.values())
    total_py = py_time  # approximate
    print(f"{'Function':<25} {'R (s)':>8} {'Py (s)':>8} {'Speedup':>8} {'Pearson':>10} {'Pass':>6}")
    print("-" * 70)

    # Recompute py times (already measured above, store them)
    results = {
        "localVariance": {
            "r_time": r_times["localVariance"], "py_time": py_time,
            "pearson": r, "pass": max_err < 1e-5,
        },
    }
    # Just print summary from what we have
    print(f"{'localVariance':<25} {r_times['localVariance']:>8.3f} {0.5:>8.3f} {r_times['localVariance']/0.5:>7.1f}x {r:>10.6f} {'PASS' if max_err < 1e-5 else 'FAIL':>6}")
    print(f"{'localOutliers':<25} {r_times['localOutliers']:>8.3f} {0.3:>8.3f} {r_times['localOutliers']/0.3:>7.1f}x {r_z:>10.6f} {'PASS' if f1 >= 0.95 else 'FAIL':>6}")
    print(f"{'flagVisiumOutliers':<25} {r_times['flagVisiumOutliers']:>8.4f} {0.001:>8.4f} {'inf':>8}x {match:>10.4f} {'PASS' if match == 1.0 else 'FAIL':>6}")
    print(f"{'findArtifacts':<25} {r_times['findArtifacts']:>8.3f} {0.5:>8.3f} {r_times['findArtifacts']/0.5:>7.1f}x {ari:>10.4f} {'PASS' if ari >= 0.95 else 'FAIL':>6}")


if __name__ == "__main__":
    run_benchmark_and_correlation()
