import numpy as np
import pandas as pd
from omicsqc.models import ProteinData

def filter_data(
    protdata: ProteinData,
    protein_missingness_threshold_pct: np.float16,
    sample_missingness_threshold_pct: np.float16,
    batch_col: str | None = None,
):
    missingness = protdata.missingness_report
    valid_proteins = missingness.per_protein_missingness[
        missingness.per_protein_missingness <= protein_missingness_threshold_pct
    ].index
    valid_samples = missingness.per_sample_missingness[
        missingness.per_sample_missingness <= sample_missingness_threshold_pct
    ].index

    filtered_data = protdata.data.loc[valid_proteins, valid_samples]
    filtered_metadata = protdata.metadata.loc[valid_samples]

    if batch_col is not None:
        fully_missing_in_any_batch = pd.Series(False, index=valid_proteins)
        for batch in filtered_metadata[batch_col].dropna().unique():
            batch_samples = filtered_metadata.index[
                filtered_metadata[batch_col] == batch
            ]
            batch_samples = batch_samples.intersection(filtered_data.columns)
            if len(batch_samples) == 0:
                continue
            all_nan = filtered_data.loc[:, batch_samples].isna().all(axis=1)
            fully_missing_in_any_batch |= all_nan

        n_dropped = fully_missing_in_any_batch.sum()
        if n_dropped > 0:
            protdata.logs.append(
                f"Dropped {n_dropped} proteins fully missing within at least "
                f"one level of '{batch_col}' (global missingness passed "
                f"threshold but per-batch completeness failed)."
            )
        valid_proteins = valid_proteins[~fully_missing_in_any_batch]
        filtered_data = filtered_data.loc[valid_proteins]
    protdata.invalid_proteins = fully_missing_in_any_batch.index.to_list()
    protdata.valid_proteins = valid_proteins
    protdata.filtered_data = filtered_data
    protdata.filtered_metadata = filtered_metadata
    return protdata

def median_centering(protdata:ProteinData):
    
    log2_expr = np.log2(protdata.filtered_data)
    
    sample_medians = log2_expr.median(axis=0)  
    normalized = log2_expr - sample_medians 
    
    protdata.transformed_data = normalized
    protdata.logs.append("Log2 transform and median normalization applied")
    return protdata


def transform_data(
    protdata: ProteinData,
    protein_missingness_threshold_pct,
    sample_missingness_threshold_pct,
    batch_col: str | None = None,
):
    protdata = filter_data(
        protdata,
        protein_missingness_threshold_pct,
        sample_missingness_threshold_pct,
        batch_col=batch_col,
    )
    protdata = median_centering(protdata)
    return protdata

