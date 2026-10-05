import logging
import warnings

import numpy as np
import pandas as pd
from scipy.stats import f as f_dist, false_discovery_control, kruskal

from ..config import QCconfigs
from ..models import BatchMetrics

logger = logging.getLogger(__name__)


def _aligned_batch(data_df: pd.DataFrame, batch: pd.Series):
    """Align batch labels to the expression columns and drop samples with no label.
    Returns (values as float array [proteins x kept samples], labels array)."""
    batch = batch.reindex(data_df.columns)
    keep = batch.notna().to_numpy()
    values = data_df.to_numpy(dtype=float)[:, keep]
    labels = batch.to_numpy()[keep]
    return values, labels

# Benjamini-Hochberg for controlling False Discovery rate (less strict than Bonferri correction)
def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values (untestable proteins) are skipped."""
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    if ok.any():
        q[ok] = false_discovery_control(p[ok], method="bh")
    return q


# ANOVA test Vectorized
def vectorized_anova(data_df: pd.DataFrame, batch: pd.Series,
                     alpha: float = 0.05, ANOVA_eta_threshold: float = 0.1,
                     min_batches: int = 3) -> pd.DataFrame:
    
    if min_batches < 2:
        raise ValueError("min_batches must be at least 2 (ANOVA needs two groups).")
    expr, labels = _aligned_batch(data_df, batch)
 
    if np.nanmax(expr) > 100:
        raise ValueError("Max value > 100: data does not look log2-transformed.")
    if np.nanmax(np.abs(np.nanmedian(expr, axis=0))) > 0.1:
        warnings.warn("Sample medians are not near 0: data may not be median-centered.")
 
    n_prot = expr.shape[0]
    nan = lambda: np.full(n_prot, np.nan)
 
    n = (~np.isnan(expr)).sum(axis=1)
    grand_mean = np.divide(np.nansum(expr, axis=1), n, out=nan(), where=n > 0)
 
    ss_between = np.zeros(n_prot)
    ss_within = np.zeros(n_prot)
    k = np.zeros(n_prot, dtype=int)                
 
    for b in np.unique(labels):
        group = expr[:, labels == b]
        gn = (~np.isnan(group)).sum(axis=1)
        has = gn > 0
        k += has
        gmean = np.divide(np.nansum(group, axis=1), gn, out=nan(), where=has)
        ss_between[has] += gn[has] * (gmean[has] - grand_mean[has]) ** 2
        ss_within += np.nansum((group - gmean[:, None]) ** 2, axis=1)
 
    df_b, df_w = k - 1, n - k
    statistically_ok = (df_b > 0) & (df_w > 0) & (ss_within > 0)
    low_coverage = k < min_batches                   
    ok = statistically_ok & ~low_coverage             
 
    excluded = data_df.index[statistically_ok & low_coverage]
    if len(excluded):
        logger.info("ANOVA: excluded %d proteins seen in fewer than %d batches "
                    "(first 10: %s). Full list: result[result['low_coverage']].",
                    len(excluded), min_batches, list(excluded[:10]))
        logger.debug("ANOVA excluded proteins: %s", list(excluded))
 
    F, p, eta = nan(), nan(), nan()
    F[ok] = (ss_between[ok] / df_b[ok]) / (ss_within[ok] / df_w[ok])
    p[ok] = f_dist.sf(F[ok], df_b[ok], df_w[ok])
    eta[ok] = ss_between[ok] / (ss_between[ok] + ss_within[ok])
    q = _bh(p)
 
    out = pd.DataFrame({"F": F, "p_value": p, "q_value": q, "eta_sq": eta,
                        "n_batches": k, "low_coverage": low_coverage}, index=data_df.index)
    out["significant"] = out["q_value"] < alpha
    out["Flag"] = out["significant"] & (out["eta_sq"] >= ANOVA_eta_threshold)
    logger.info("ANOVA executed on %d proteins (%d testable, min_batches=%d)",
                n_prot, int(ok.sum()), min_batches)
    return out.sort_values("F", ascending=False)

# Kruskal-Wallis test
def kruskal_batch(data_df: pd.DataFrame, batch: pd.Series,
                  alpha: float = QCconfigs.alpha) -> pd.DataFrame:

    values, labels = _aligned_batch(data_df, batch)
    masks = [labels == b for b in np.unique(labels)]

    n_prot = values.shape[0]
    H = np.full(n_prot, np.nan)
    P = np.full(n_prot, np.nan)

    for i in range(n_prot):
        protein = values[i]
        groups = []
        for mask in masks:
            g = protein[mask]
            g = g[~np.isnan(g)]
            if len(g):
                groups.append(g)
        if len(groups) < 2:
            continue
        try:
            H[i], P[i] = kruskal(*groups)
        except ValueError:
            pass

    result = pd.DataFrame({"H": H, "p_value": P, "q_value": _bh(P)}, index=data_df.index)
    result["significant"] = result["q_value"] < alpha
    logger.info("Kruskal-Wallis executed on %d proteins (%d testable)",
                n_prot, int((~np.isnan(P)).sum()))
    return result


def batch_variation_test(data_df: pd.DataFrame, meta_df: pd.DataFrame, batch_col: str,
                         batch_effect_method: str = QCconfigs.batch_effect_method,
                         alpha: float = QCconfigs.alpha, eta_threshold: float = 0.1):
    """Run the per-protein batch test.
    data_df must be log2-transformed and median-centered (e.g. protdata.transformed_data
    or protdata.corrected_data), never the raw protdata.data."""
    batch = meta_df[batch_col]

    if batch_effect_method == "anova":
        per_protein = vectorized_anova(data_df, batch, alpha=alpha,
                                       ANOVA_eta_threshold=eta_threshold)
    elif batch_effect_method == "kruskal_wallis":
        per_protein = kruskal_batch(data_df, batch, alpha=alpha)
    else:
        raise ValueError(f"Unknown batch_effect_method '{batch_effect_method}'. "
                         "Use 'anova' or 'kruskal_wallis'.")

    n_testable = int(per_protein["p_value"].notna().sum())
    n_sig = int(per_protein["significant"].sum())
    n_flag = int(per_protein["Flag"].sum()) if "Flag" in per_protein else None

    def pct(x):
        return 100 * x / n_testable if (x is not None and n_testable) else np.nan

    batch_metrics = BatchMetrics(
        batch_effect_per_protein=per_protein,
        n_testable=n_testable,
        n_significant=n_sig,
        n_significant_pct=pct(n_sig),
        n_flagged=n_flag,
        n_flagged_pct=pct(n_flag),
    )

    lines = [
        "[Batch Variation Test]",
        f"Statistical method:          {batch_effect_method}",
        f"Significance level (alpha):  {alpha} (BH-adjusted)",
        f"Testable proteins:           {n_testable}",
        f"Significant proteins:        {n_sig} ({pct(n_sig):.2f}%)",
    ]
    if n_flag is not None:
        lines.append(f"Flagged (q<alpha, eta^2>={eta_threshold}): {n_flag} ({pct(n_flag):.2f}%)")
    logs = "\n".join(lines)

    logger.info("Batch variation test (%s): %d/%d significant", batch_effect_method, n_sig, n_testable)
    return batch_metrics, logs