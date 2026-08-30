# config.py
import numpy as np
from dataclasses import dataclass
@dataclass

class QCconfigs:
    missingness_pval_threshold: np.float32 = 0.05
    missingness_corr_threshold: np.float32 = 0.3
    batch_effect_method: str = "anova"
    alpha: np.float32 = 0.05
    MinProb_lmbda: float = 1.8,
    MinProb_delta: float = 0.3,
