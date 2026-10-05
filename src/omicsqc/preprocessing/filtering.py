import logging
import pandas as pd
logger = logging.getLogger(__name__)


def filter_data(df: pd.DataFrame,
                max_protein_missing_pct: float = 80.0,
                max_sample_missing_pct: float = 50.0,
                batch: pd.Series | None = None):
    """Drop samples, then proteins, with too much missingness.

    Lenient by default on purpose: low-abundance proteins are the ones with most missing
    values, and they are exactly what the imputation step is for. A strict filter leaves
    almost nothing to impute.

    Samples are filtered first, and protein missingness is computed on the samples that remain.
    Proteins with no observed value at all are always dropped.
    batch (optional, labels indexed by sample id) is only used to warn when dropping samples
    leaves a batch with fewer than 2 samples or removes it entirely.

    Returns (filtered_df, protein_mask, sample_mask, logs). The masks are boolean Series over
    the input's proteins and samples (True = kept)."""
    
    for name, value in [("max_protein_missing_pct", max_protein_missing_pct),
                        ("max_sample_missing_pct", max_sample_missing_pct)]:
        if not 0 <= value <= 100:
            raise ValueError(f"{name} must be between 0 and 100, got {value}.")

    sample_missing = df.isna().mean(axis=0) * 100
    sample_mask = sample_missing <= max_sample_missing_pct
    kept = df.loc[:, sample_mask]

    protein_missing = kept.isna().mean(axis=1) * 100
    protein_mask = (protein_missing <= max_protein_missing_pct) & kept.notna().any(axis=1)
    filtered = kept.loc[protein_mask]

    if batch is not None:
        before = batch.reindex(df.columns).dropna()
        after = batch.reindex(filtered.columns).dropna().value_counts()
        lost = sorted(set(before) - set(after.index))
        small = after[after < 2].index.tolist()
        if lost:
            logger.warning("Sample filtering removed every sample of batches: %s", lost)
        if small:
            logger.warning("Batches left with fewer than 2 samples: %s", small)

    pct_before = df.isna().mean().mean() * 100
    pct_after = filtered.isna().mean().mean() * 100 if filtered.size else float("nan")
    logs = "\n".join([
        "[Filtering]",
        f"Thresholds: proteins <= {max_protein_missing_pct}% missing, samples <= {max_sample_missing_pct}% missing",
        f"Samples kept:   {int(sample_mask.sum())} of {len(sample_mask)}",
        f"Proteins kept:  {int(protein_mask.sum())} of {len(protein_mask)}",
        f"Missing values: {pct_before:.2f}% before, {pct_after:.2f}% after",
    ])
    logger.info("Filtering kept %d proteins and %d samples", int(protein_mask.sum()), int(sample_mask.sum()))
    return filtered, protein_mask, sample_mask, logs