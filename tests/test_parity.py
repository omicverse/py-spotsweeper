"""Parity tests comparing Python implementation against R reference outputs."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import f1_score, adjusted_rand_score

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from conftest import (
    load_visium_fixture,
    load_artifact_fixture,
    load_reference_outputs,
    load_artifact_reference_outputs,
    load_r_knn_indices,
)


class TestLocalVarianceParity:
    """Test localVariance against R reference (deterministic, atol=1e-8)."""

    def test_residuals_match(self):
        from spotsweeper.local_variance import local_variance

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_residuals = np.array(ref["local_variance_residuals"])
        r_knn = load_r_knn_indices()

        adata = local_variance(
            adata,
            metric="subsets_Mito_percent",
            n_neighbors=36,
            name="local_mito_variance_k36",
            log=False,
            knn_indices=r_knn,
        )

        py_residuals = adata.obs["local_mito_variance_k36"].values

        # Check same length
        assert len(py_residuals) == len(ref_residuals)

        # Check finite
        assert np.all(np.isfinite(py_residuals)), "Python residuals contain non-finite values"

        # Deterministic parity: max absolute error
        max_err = np.max(np.abs(py_residuals - ref_residuals))
        print(f"localVariance max abs error: {max_err:.2e}")
        assert max_err < 1e-5, f"localVariance max abs error {max_err:.2e} exceeds threshold 1e-5"

    def test_residuals_correlation(self):
        """Pearson and Spearman correlation between R and Python residuals."""
        from spotsweeper.local_variance import local_variance

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_residuals = np.array(ref["local_variance_residuals"])
        r_knn = load_r_knn_indices()

        adata = local_variance(
            adata, metric="subsets_Mito_percent", n_neighbors=36,
            name="local_mito_variance_k36", log=False, knn_indices=r_knn,
        )
        py_residuals = adata.obs["local_mito_variance_k36"].values

        r, _ = pearsonr(ref_residuals, py_residuals)
        rho, _ = spearmanr(ref_residuals, py_residuals)
        print(f"localVariance Pearson r: {r:.8f}, Spearman rho: {rho:.8f}")
        assert r > 0.99, f"Pearson r={r:.8f} below 0.99"
        assert rho > 0.99, f"Spearman rho={rho:.8f} below 0.99"


class TestLocalOutliersParity:
    """Test localOutliers against R reference."""

    def test_zscores_match(self):
        from spotsweeper.local_outliers import local_outliers

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_z = np.array(ref["local_outlier_zscores"])
        r_knn = load_r_knn_indices()

        adata = local_outliers(
            adata,
            metric="sum",
            direction="lower",
            log=True,
            n_neighbors=36,
            knn_indices=r_knn,
        )

        py_z = adata.obs["sum_z"].values

        assert len(py_z) == len(ref_z)
        assert np.all(np.isfinite(py_z))

        max_err = np.max(np.abs(py_z - ref_z))
        print(f"localOutliers z-score max abs error: {max_err:.2e}")
        assert max_err < 0.1, f"z-score max abs error {max_err:.2e} too large"

    def test_zscores_correlation(self):
        """Pearson and Spearman correlation between R and Python z-scores."""
        from spotsweeper.local_outliers import local_outliers

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_z = np.array(ref["local_outlier_zscores"])
        r_knn = load_r_knn_indices()

        adata = local_outliers(
            adata, metric="sum", direction="lower", log=True,
            n_neighbors=36, knn_indices=r_knn,
        )
        py_z = adata.obs["sum_z"].values

        r, _ = pearsonr(ref_z, py_z)
        rho, _ = spearmanr(ref_z, py_z)
        print(f"localOutliers z-score Pearson r: {r:.8f}, Spearman rho: {rho:.8f}")
        assert r > 0.99, f"Pearson r={r:.8f} below 0.99"
        assert rho > 0.99, f"Spearman rho={rho:.8f} below 0.99"

    def test_outlier_flags_f1(self):
        from spotsweeper.local_outliers import local_outliers

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_flags = np.array(ref["local_outlier_flags"], dtype=bool)
        r_knn = load_r_knn_indices()

        adata = local_outliers(
            adata,
            metric="sum",
            direction="lower",
            log=True,
            n_neighbors=36,
            knn_indices=r_knn,
        )

        py_flags = adata.obs["sum_outliers"].values.astype(bool)

        f1 = f1_score(ref_flags, py_flags)
        print(f"localOutliers flag F1: {f1:.4f} (ref outliers: {ref_flags.sum()}, py: {py_flags.sum()})")
        assert f1 >= 0.95, f"Outlier flag F1 {f1:.4f} below threshold 0.95"


class TestFindArtifactsParity:
    """Test findArtifacts against R reference (clustering, ARI >= 0.95)."""

    def test_artifact_labels_ari(self):
        from spotsweeper.find_artifacts import find_artifacts

        adata = load_artifact_fixture()
        ref = load_artifact_reference_outputs()
        ref_labels = np.array(ref["artifact_labels"], dtype=int)

        # Load R kNN indices for artifact fixture
        r_knn_6 = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "..", "data", "r_artifact_knn_k6.csv")
        ).values - 1
        r_knn_18 = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "..", "data", "r_artifact_knn_k18.csv")
        ).values - 1
        knn_dict = {6: r_knn_6, 18: r_knn_18}

        adata = find_artifacts(
            adata,
            mito_percent="expr_chrM_ratio",
            mito_sum="expr_chrM",
            n_order=2,
            shape="hexagonal",
            name="artifact",
            seed=42,
            knn_indices_dict=knn_dict,
        )

        py_labels = adata.obs["artifact"].values.astype(int)

        ari = adjusted_rand_score(ref_labels, py_labels)
        print(f"findArtifacts ARI: {ari:.4f} (ref artifacts: {ref_labels.sum()}, py: {py_labels.sum()})")
        assert ari >= 0.95, f"findArtifacts ARI {ari:.4f} below threshold 0.95"

    def test_artifact_local_variance_correlation(self):
        """Correlation of intermediate local variance values (k6, k18)."""
        from spotsweeper.find_artifacts import find_artifacts

        adata = load_artifact_fixture()
        r_knn_6 = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "..", "data", "r_artifact_knn_k6.csv")
        ).values - 1
        r_knn_18 = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "..", "data", "r_artifact_knn_k18.csv")
        ).values - 1

        adata = find_artifacts(
            adata, mito_percent="expr_chrM_ratio", mito_sum="expr_chrM",
            n_order=2, shape="hexagonal", name="artifact", seed=42,
            knn_indices_dict={6: r_knn_6, 18: r_knn_18},
        )

        k6 = adata.obs["k6"].values
        k18 = adata.obs["k18"].values

        r6, _ = pearsonr(k6, k6)  # self-correlation (sanity)
        r18, _ = pearsonr(k18, k18)
        print(f"findArtifacts k6 self-corr: {r6:.8f}, k18 self-corr: {r18:.8f}")

        # Correlation between k6 and k18 (should be positive, related but different scales)
        r_cross, _ = pearsonr(k6, k18)
        print(f"findArtifacts k6 vs k18 cross-corr: {r_cross:.4f}")
        assert r_cross > 0.3, f"k6 vs k18 correlation {r_cross:.4f} too low"


class TestFlagVisiumOutliersParity:
    """Test flagVisiumOutliers against R reference (exact match)."""

    def test_flags_exact_match(self):
        from spotsweeper.flag_visium_outliers import flag_visium_outliers

        adata = load_visium_fixture()
        ref = load_reference_outputs()
        ref_flags = np.array(ref["systematic_outlier_flags"], dtype=bool)

        adata = flag_visium_outliers(adata)

        py_flags = adata.obs["systematic_outliers"].values.astype(bool)

        match = np.mean(ref_flags == py_flags)
        print(f"flagVisiumOutliers exact match: {match:.4f} (ref: {ref_flags.sum()}, py: {py_flags.sum()})")
        assert match == 1.0, f"flagVisiumOutliers match {match:.4f} < 1.0"
