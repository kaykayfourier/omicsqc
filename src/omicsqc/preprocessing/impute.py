import numpy as np
import pandas as pd


def MinProbImpute(
    data_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    batch_col: str,
    lmbda: float = 1.8,
    delta: float = 0.3,
    random_state: int | None = None,
) -> pd.DataFrame:
    
    rng = np.random.default_rng(random_state)
    imputed = data_df.astype(np.float64).copy()

    full_matrix = data_df.to_numpy(dtype=float)
    global_mu = np.nanmean(full_matrix, axis=1, keepdims=True)
    global_sigma = np.nanstd(full_matrix, axis=1, keepdims=True, ddof=1)

    for batch in meta_df[batch_col].dropna().unique():
        batch_samples = meta_df.index[meta_df[batch_col] == batch]
        batch_samples = batch_samples.intersection(imputed.columns)
        if len(batch_samples) == 0:
            continue

        matrix = imputed.loc[:, batch_samples].to_numpy(dtype=float, copy = True)
        mask = np.isnan(matrix)
        if not mask.any():
            continue

        mu = np.nanmean(matrix, axis=1, keepdims=True)
        sigma = np.nanstd(matrix, axis=1, keepdims=True, ddof=1)

        # fallback: protein fully missing in this batch -> use global stats
        empty_rows = np.isnan(mu).flatten()
        if empty_rows.any():
            mu[empty_rows] = global_mu[empty_rows]
            sigma[empty_rows] = global_sigma[empty_rows]

        # guard degenerate sigma (NaN, 0, or single-observation batches)
        bad_sigma = np.isnan(sigma) | (sigma <= 0)
        if bad_sigma.any():
            valid_sigma = sigma[~bad_sigma]
            fallback_val = (
                np.median(valid_sigma) if valid_sigma.size > 0 else 1e-6
            )
            sigma[bad_sigma] = fallback_val

        mu_imp = mu - lmbda * sigma
        sigma_imp = delta * sigma

        random_matrix = rng.normal(
            loc=mu_imp, scale=sigma_imp, size=matrix.shape
        )
        matrix[mask] = random_matrix[mask]
        imputed.loc[:, batch_samples] = matrix

    return imputed