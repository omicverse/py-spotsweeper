"""SpotSweeper: Spatially-aware quality control for spatial transcriptomics."""

from spotsweeper.local_variance import local_variance
from spotsweeper.local_outliers import local_outliers
from spotsweeper.find_artifacts import find_artifacts
from spotsweeper.flag_visium_outliers import flag_visium_outliers

__version__ = "0.1.0"
