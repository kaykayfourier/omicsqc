import logging
import warnings

import numpy as np
import pandas as pd

from ..models import VarianceReport

logger = logging.getLogger(__name__)


# Adjusted R^2 of one grouping factor for every protein: how much of the protein's variance the groups explain, corrected for the number of groups.
# Plain R^2 grows with the number of groups even for pure noise (about (k-1)/(n-1)), so a
# 67-level batch would look bigger than a 2-level diagnosis. The adjustment makes them comparable.
# data_df must be log2-transformed and median-centered. groups is a Series indexed by sample id.
def variance_explained(data_df: pd.DataFrame, groups: pd.Series, min_groups: int = 2, min_obs: int = 30):
    factor = groups.name
    if min_groups < 2:
        raise ValueError("min_groups must be at least 2.")
    if pd.api.types.is_numeric_dtype(groups) and groups.nunique() > 10:
        raise ValueError(f"Factor '{factor}' looks continuous ({groups.nunique()} distinct numeric values). "
                         "Only categorical factors are supported.")

    if data_df.max().max() > 100:
        raise ValueError("Max value > 100: data does not look log2-transformed.")
    if data_df.median(axis=0).abs().max() > 0.1:
        warnings.warn("Sample medians are not near 0: data may not be median-centered.")

    groups = groups.reindex(data_df.columns)
    keep = groups.notna().to_numpy()
    expr = data_df.to_numpy(dtype=float)[:, keep]
    labels = groups.to_numpy()[keep]
    if len(np.unique(labels)) < 2:
        raise ValueError(f"Factor '{factor}' has fewer than 2 groups among the samples in the data.")

    n_prot = expr.shape[0]
    nan = lambda: np.full(n_prot, np.nan)

    n = (~np.isnan(expr)).sum(axis=1)
    grand_mean = np.divide(np.nansum(expr, axis=1), n, out=nan(), where=n > 0)

    ss_between = np.zeros(n_prot)
    ss_within = np.zeros(n_prot)
    k = np.zeros(n_prot, dtype=int)                 

    for g in np.unique(labels):
        group = expr[:, labels == g]
        gn = (~np.isnan(group)).sum(axis=1)
        has = gn > 0
        k += has
        gmean = np.divide(np.nansum(group, axis=1), gn, out=nan(), where=has)
        ss_between[has] += gn[has] * (gmean[has] - grand_mean[has]) ** 2
        ss_within += np.nansum((group - gmean[:, None]) ** 2, axis=1)

    ss_total = ss_between + ss_within
    
    ok = (k >= min_groups) & (n >= min_obs) & (n - k > 0) & (ss_total > 0)

    # adjusted R^2 = 1 - (ss_within / (n - k)) / (ss_total / (n - 1)). Can be negative:
    # that means the groups explain no more than chance.
    adj_r2 = nan()
    adj_r2[ok] = 1 - (ss_within[ok] / (n[ok] - k[ok])) / (ss_total[ok] / (n[ok] - 1))

    n_excluded = int((~ok).sum())
    if n_excluded:
        logger.info("Variance explained by '%s': %d of %d proteins not testable "
                    "(fewer than %d groups or %d samples observed, or no variance).",
                    factor, n_excluded, n_prot, min_groups, min_obs)

    variance_per_protein = pd.Series(adj_r2, index=data_df.index, name=factor)
    summary = {
        "factor": factor,
        "n_groups": len(np.unique(labels)),
        "n_samples": int(keep.sum()),                  # samples this factor was computed on
        "n_proteins_tested": int(ok.sum()),
        "mean_r2": variance_per_protein.mean(),        # these are adjusted R^2 values
        "median_r2": variance_per_protein.median(),
        "std_r2": variance_per_protein.std(),
        "q25": variance_per_protein.quantile(0.25),
        "q75": variance_per_protein.quantile(0.75),
        "min_r2": variance_per_protein.min(),
        "max_r2": variance_per_protein.max(),
    }
    return variance_per_protein, summary


# Runs the same calculation for several factors so batch can be compared with the biology.
# complete_cases=True: every factor is computed on the same samples, the ones that have a value for all factors. 
# Otherwise a factor that is only recorded for part of the data (e.g. a clinical score missing for one cohort) is measured on a different set of samples than the others and
# the numbers are not comparable. The cost is dropping samples, so the count is logged.
def compare_explained_variances(data_df: pd.DataFrame, meta_df: pd.DataFrame, factors: list[str],
                                min_groups: int = 2, min_obs: int = 30, complete_cases: bool = True):
    for factor in factors:
        if factor not in meta_df.columns:
            raise ValueError(f"Factor '{factor}' is not a column of the metadata.")

    n_before = data_df.shape[1]
    if complete_cases:
        shared = data_df.columns.intersection(meta_df.index)
        has_all = meta_df.loc[shared, factors].notna().all(axis=1)
        keep = has_all.index[has_all]
        if len(keep) == 0:
            raise ValueError("No samples have a value for every factor.")
        data_df = data_df[keep]
        logger.info("Complete cases: using %d of %d samples.", len(keep), n_before)

    summaries = []
    per_protein_frames = []
    for factor in factors:
        r2_series, summary = variance_explained(data_df, meta_df[factor],
                                                min_groups=min_groups, min_obs=min_obs)
        summaries.append(summary)
        per_protein_frames.append(r2_series)

    comparison_df = pd.concat(per_protein_frames, axis=1)
    summary_df = pd.DataFrame(summaries).set_index("factor")

    variation_test_report = VarianceReport(factors=factors, comparison_df=comparison_df, summary_df=summary_df)

    logs = "\n".join([
        "[Explained Variance Comparison] (adjusted R^2)",
        f"Factors analyzed:             {', '.join(factors)}",
        f"Per-protein comparison shape: {comparison_df.shape}",
        f"Samples used:                {data_df.shape[1]} of {n_before}",
        summary_df[["n_groups", "n_samples", "n_proteins_tested", "median_r2"]].to_string(),
    ])
    logger.info("Explained variance computed for %d factors", len(factors))
    return variation_test_report, logs