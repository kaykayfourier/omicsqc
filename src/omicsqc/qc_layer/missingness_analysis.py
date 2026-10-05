import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kruskal

from ..config import QCconfigs
from ..models import MissingnessReport

logger = logging.getLogger(__name__)


# In TMT a protein is usually either identified in a whole plex or not at all, so
# detection is counted per plex (batch), not per sample.
def detection_per_plex(data_df, batch):
    seen = data_df.notna().T.groupby(batch.values).any().T   # protein x plex, True if seen in any sample
    return seen.mean(axis=1)                                 # fraction of plexes the protein was seen in


# Does detection go up with protein intensity? (a sign of left-censoring)
# We can only say yes/no here. MAR vs MNAR can't be told apart from observed data.
def intensity_dependence(log2_df, batch, corr_threshold=QCconfigs.missingness_corr_threshold):
    detection = detection_per_plex(log2_df, batch)
    mean_intensity = log2_df.mean(axis=1)

    keep = detection > 0                      # proteins never seen have no mean intensity
    detection = detection[keep]
    mean_intensity = mean_intensity[keep]

    rho, pval = spearmanr(detection, mean_intensity)

    # detection rate per intensity decile, easier to read than a single rho
    bins = pd.qcut(mean_intensity, 10, duplicates="drop")
    curve = detection.groupby(bins, observed=True).mean()
    curve.index = [f"{iv.left:.1f} to {iv.right:.1f}" for iv in curve.index]

    if np.isnan(rho):
        verdict = "Intensity dependence could not be computed (detection rate is constant)."
    elif rho > corr_threshold:
        verdict = "Intensity-dependent missingness: yes (detection rises with intensity, consistent with left-censoring)."
    elif rho < -corr_threshold:
        verdict = "Intensity-dependent missingness: unexpected negative association, check the data."
    else:
        verdict = "Intensity-dependent missingness: no strong association (MCAR not rejected on this axis)."
    verdict += " MAR vs MNAR cannot be distinguished from observed data."

    return rho, pval, curve, verdict, len(detection)


# Share of (protein, plex) cells that are only partly missing.
# Near 0 means proteins are missing for whole plexes (identification level).
def plex_structure(data_df, batch):
    frac = data_df.notna().T.groupby(batch.values).mean().T   # fraction of samples detected, per protein per plex
    present = frac > 0
    partial = present & (frac < 1)
    return partial.sum().sum() / present.sum().sum()


# Does sample missingness depend on an observed covariate (batch, cohort, age...)?
# Outcome is each sample's missingness %. Categorical -> Kruskal-Wallis, numeric -> Spearman.
def covariate_missingness(sample_missing, meta_df, covariates):
    for col in covariates:
        if col not in meta_df.columns:
            raise ValueError(f"Covariate '{col}' is not a column of the metadata.")

    rows = []
    for col in covariates:
        values = meta_df[col].reindex(sample_missing.index)
        ok = values.notna()

        if ok.sum() < len(values) / 2:
            logger.warning("Covariate '%s' skipped: more than half the values are missing.", col)
            continue
        values = values[ok]
        miss = sample_missing[ok]

        if values.nunique() < 2:
            logger.warning("Covariate '%s' skipped: only one level.", col)
            continue

        # numeric with many distinct values is treated as continuous, everything else as groups
        if pd.api.types.is_numeric_dtype(values) and values.nunique() > 10:
            test = "spearman"
            stat, pval = spearmanr(values, miss)
            effect = stat
        else:
            test = "kruskal"
            groups = [miss[values == v].to_numpy() for v in values.unique()]
            n, k = len(values), len(groups)
            if n - k < 1:
                logger.warning("Covariate '%s' skipped: too many levels for its sample count.", col)
                continue
            try:
                stat, pval = kruskal(*groups)
            except ValueError:
                logger.warning("Covariate '%s' skipped: Kruskal-Wallis could not be computed.", col)
                continue
            effect = max((stat - k + 1) / (n - k), 0)    # eta-squared based on H

        rows.append({"covariate": col, "test": test, "statistic": stat,
                     "effect_size": effect, "p_value": pval, "n_samples": len(values)})

    return pd.DataFrame(rows, columns=["covariate", "test", "statistic",
                                       "effect_size", "p_value", "n_samples"])


# data_df: raw intensities (not log2). Zeros are treated as missing.
# covariates: metadata columns to test against missingness. Defaults to just the batch column.
def compute_missingness(data_df, meta_df, batch_col, covariates=None,
                        corr_threshold=QCconfigs.missingness_corr_threshold):

    if batch_col not in meta_df.columns:
        raise ValueError(f"Given batch_col: {batch_col} does not exist in metadata columns.")
    if covariates is None:
        covariates = [batch_col]

    data_df = data_df.replace(0.0, np.nan)

    overall_missingness = data_df.isna().mean().mean() * 100
    if overall_missingness < 1e-8:
        return None, "missingness analysis not applicable because no missing values remain"

    per_protein_missingness = data_df.isna().mean(axis=1) * 100
    per_sample_missingness = data_df.isna().mean(axis=0) * 100

    # only samples that have a batch label are used for the batch based numbers
    batch = meta_df[batch_col].reindex(data_df.columns)
    labeled = data_df.loc[:, batch.notna()]
    batch = batch[batch.notna()]

    batch_missingness = (labeled.isna().T.groupby(batch.values).mean().mean(axis=1) * 100
                         ).sort_values(ascending=False)

    log2_df = np.log2(labeled)
    rho, pval, curve, verdict, n_tested = intensity_dependence(log2_df, batch, corr_threshold)
    partial_fraction = plex_structure(labeled, batch)
    covariate_table = covariate_missingness(per_sample_missingness, meta_df, covariates)

    report = MissingnessReport(
        per_protein_missingness=per_protein_missingness,
        per_sample_missingness=per_sample_missingness,
        overall_missingness=overall_missingness,
        batch_missingness=batch_missingness,
        spearman_r=rho,
        spearman_pval=pval,
        missingness_inference=verdict,
        detection_curve=curve,
        n_proteins_tested=n_tested,
        partial_fraction=partial_fraction,
        covariate_table=covariate_table,
    )

    lines = [
        "[Missingness Analysis]",
        f"Per-protein missingness (mean): {per_protein_missingness.mean():.2f}%",
        f"Per-sample missingness (mean):  {per_sample_missingness.mean():.2f}%",
        f"Overall missingness:            {overall_missingness:.2f}%",
        f"Spearman (detection vs intensity): r = {rho:.4f}, p = {pval:.2e}, n = {n_tested} proteins",
        f"Partly-missing protein-plex cells: {partial_fraction * 100:.1f}%",
        verdict,
    ]
    if len(covariate_table):
        lines.append("Covariate vs sample missingness:")
        lines.append(covariate_table.to_string(index=False))
    logs = "\n".join(lines)

    logger.info("Missingness analysis done (%.2f%% missing overall)", overall_missingness)
    return report, logs