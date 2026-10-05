import warnings
import numpy as np
import pandas as pd


def log2_median_centering(df: pd.DataFrame) -> pd.DataFrame:
    """log2-transform raw intensities, then median-center each sample (column).
    Zeros are treated as missing and NaNs are preserved.
    Refuses data that already looks log2-transformed, so it can't be applied twice."""
    values = df.to_numpy(dtype=float)
    if np.isnan(values).all():
        raise ValueError("The matrix is entirely missing.")

    # raw intensities can't be negative, but centered log2 values usually are
    if np.nanmin(values) < 0:
        raise ValueError("Negative values found: the data looks already log2-transformed and centered.")

    # raw TMT intensities run into the thousands or more; log2 values stay below ~40
    if np.nanmax(values) < 50:
        raise ValueError(f"Max value is {np.nanmax(values):.1f}: the data looks already log2-transformed.")

    n_zero = int((values == 0).sum())
    if n_zero:
        warnings.warn(f"{n_zero} zero values treated as missing before log2.")

    log2 = np.log2(df.astype(float).replace(0.0, np.nan))
    return log2 - log2.median(axis=0)       # pandas median skips NaNs