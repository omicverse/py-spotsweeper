# MATH.md — Perturbation bounds

## Class A translation (no bounded ε-approximation rewrites)

py-SpotSweeper is a **Class A** port — direct translation of R algorithms without acceleration rewrites. No (B) bounded ε-approximation rewrites were applied.

The remaining numerical difference (max abs error ≈ 2.67e-6 for localVariance) comes from:

1. **IRLS convergence**: The Huber M-estimator (IRLS) converges to slightly different solutions due to different floating-point accumulation in R's MASS::rlm vs the Python implementation. This is inherent cross-language BLAS divergence, not a deliberate approximation.

2. **kNN tie-breaking**: R's BiocNeighbors::findKNN and Python's scipy.spatial.cKDTree break tied distances differently (~0.36% of spots have different neighbor sets). Parity tests use R-exported kNN indices; production code uses Python's cKDTree.

No formal perturbation bound is required for Class A ports.
