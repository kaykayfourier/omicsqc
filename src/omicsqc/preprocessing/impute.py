import logging

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

logger = logging.getLogger(__name__)


# All three functions take a filtered, log2-transformed, median-centered matrix
# (proteins x samples) and return a new matrix with no missing values.
def _check_input(df: pd.DataFrame):
    if np.nanmax(df.to_numpy(dtype=float)) > 100:
        raise ValueError("Max value > 100: data does not look log2-transformed.")
    empty_samples = df.columns[df.isna().all(axis=0)]
    if len(empty_samples):
        raise ValueError(f"{len(empty_samples)} samples have no observed values. Filter them out first.")
    empty_proteins = df.index[df.isna().all(axis=1)]
    if len(empty_proteins):
        raise ValueError(f"{len(empty_proteins)} proteins have no observed values. Filter them out first.")


# MinProb (Lazar et al., imputeLCMD): for left-censored values.
# Each sample's missing values are drawn from a normal centered at a low quantile (q) of that sample's observed values.
# The spread is the median standard deviation of the well-observed proteins (more than half the samples observed), times tune_sigma.
def impute_minprob(df: pd.DataFrame, q: float = 0.01, tune_sigma: float = 1.0,
                   seed: int | None = None) -> pd.DataFrame:
    _check_input(df)
    values = df.to_numpy(dtype=float, copy=True)
    missing = np.isnan(values)
    if not missing.any():
        return df.copy()

    sample_low = np.nanquantile(values, q, axis=0)              # one low value per sample

    well_observed = (~missing).mean(axis=1) > 0.5
    if not well_observed.any():
        raise ValueError("No protein is observed in more than half the samples; cannot estimate the spread.")
    sd = np.nanmedian(np.nanstd(values[well_observed], axis=1, ddof=1)) * tune_sigma

    rng = np.random.default_rng(seed)
    draws = rng.normal(loc=sample_low, scale=sd, size=values.shape)   # column j uses sample_low[j]
    values[missing] = draws[missing]
    return pd.DataFrame(values, index=df.index, columns=df.columns)


# kNN: for values missing at random.
# Neighbours are proteins. A missing value (protein i, sample j) is the average of sample j
# across the k proteins most similar to protein i over the samples they share.
def impute_knn(df: pd.DataFrame, k: int = 10) -> pd.DataFrame:
    _check_input(df)
    if not df.isna().any().any():
        return df.copy()
    # sklearn treats rows as observations, so passing proteins as rows makes the neighbours proteins
    values = KNNImputer(n_neighbors=k).fit_transform(df.to_numpy(dtype=float))
    return pd.DataFrame(values, index=df.index, columns=df.columns)


# Mixed: proteins observed in fewer than detection_cutoff of the samples are treated as
# left-censored (MinProb). Proteins with only occasional gaps are treated as random (kNN).
# Returns (imputed_df, groups, logs). groups labels every protein as complete,
# left_censored or random. The imputed cells are simply df.isna().
def impute_mixed(df: pd.DataFrame, detection_cutoff: float = 0.9, k: int = 10,
                 q: float = 0.01, tune_sigma: float = 1.0, seed: int | None = None):
    _check_input(df)
    if not 0 < detection_cutoff <= 1:
        raise ValueError("detection_cutoff must be in (0, 1].")

    detection = df.notna().mean(axis=1)                # fraction of samples observed
    has_missing = detection < 1
    censored = has_missing & (detection < detection_cutoff)
    random = has_missing & ~censored

    groups = pd.Series("complete", index=df.index)
    groups[censored] = "left_censored"
    groups[random] = "random"

    imputed = df.copy()
    if censored.any():
        # sample quantiles come from the whole matrix, then only the censored proteins are kept
        imputed.loc[censored] = impute_minprob(df, q, tune_sigma, seed).loc[censored]
    if random.any():
        # neighbours are taken from complete and random proteins only. Censored proteins are
        # left out as donors because their observed values are biased upward.
        imputed.loc[random] = impute_knn(df.loc[~censored], k).loc[random]

    n_missing = df.isna().sum(axis=1)
    logs = "\n".join([
        "[Imputation: mixed]",
        f"Detection cutoff: {detection_cutoff} (below -> MinProb, at or above -> kNN)",
        f"Complete proteins:        {int((groups == 'complete').sum())}",
        f"Left-censored (MinProb):  {int(censored.sum())} proteins, {int(n_missing[censored].sum())} values",
        f"Random (kNN, k={k}):       {int(random.sum())} proteins, {int(n_missing[random].sum())} values",
    ])
    logger.info("Mixed imputation: %d MinProb, %d kNN proteins", int(censored.sum()), int(random.sum()))
    return imputed, groups, logs