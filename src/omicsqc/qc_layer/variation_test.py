import numpy as np
import pandas as pd
from ..models import VarianceReport


def variance_explained(data_df:pd.DataFrame, meta_df : pd.DataFrame, factor: str,):
    
    expr = data_df.to_numpy(dtype = float)
    factors = meta_df.loc[data_df.columns, factor].to_numpy()
    unique_groups = np.unique(factors[~pd.isna(factors)])
    n_proteins, n_samples = expr.shape

    # Safe grand mean
    total_n = np.sum(~np.isnan(expr), axis=1)
    grand_mean = np.divide(np.nansum(expr, axis=1), total_n, out=np.full(n_proteins, np.nan), where=total_n > 0)

    ss_between = np.zeros(n_proteins)
    ss_within = np.zeros(n_proteins)

    for g in unique_groups:
        mask = factors == g
        group = expr[:, mask]
        group_n = np.sum(~np.isnan(group), axis=1)
        group_mean = np.divide(np.nansum(group, axis=1), group_n, out=np.full(n_proteins, np.nan), where=group_n > 0)
        valid = group_n > 0
        ss_between[valid] += group_n[valid] * (group_mean[valid] - grand_mean[valid].squeeze())**2
        ss_within[valid] += np.nansum((group[valid] - group_mean[valid, None]) ** 2, axis=1)

    total_ss = ss_between + ss_within

    r2 = np.divide(ss_between, total_ss, out=np.full(n_proteins, np.nan), where=total_ss > 0)
    
    variance_per_protein = pd.Series(r2, index = data_df.index, name = factor )
    summary = {
        "factor": factor,
        "mean_r2": np.nanmean(r2),
        "median_r2": np.nanmedian(r2),
        "std_r2": np.nanstd(r2),
        "q25": np.nanpercentile(r2, 25),
        "q75": np.nanpercentile(r2, 75),
        "min_r2": np.nanmin(r2),
        "max_r2": np.nanmax(r2),
    }
    return variance_per_protein, summary

def compare_explained_variances(data_df:pd.DataFrame, meta_df : pd.DataFrame, factors: list[str]):
    
    summaries = []
    per_protein_frames = []

    for factor in factors:
        valid_samples = meta_df[factor].dropna().index
        valid_samples = valid_samples.intersection(data_df.columns)

        r2_series, summary = variance_explained(data_df, meta_df, factor = factor)

        summaries.append(summary)
        per_protein_frames.append(r2_series)
    comparison_df = pd.concat(per_protein_frames, axis = 1)
    summary_df = pd.DataFrame(summaries).set_index("factor")

    variation_test_report = VarianceReport(factors= factors, comparison_df= comparison_df, summary_df= summary_df)


    logs = (
        f"""
    [Explained Variance Comparison]
    Factors analyzed:               {", ".join(factors)}
    Per-protein comparison shape:   {comparison_df.shape}
    Summary dataframe shape:        {summary_df.shape}
    """.strip()
    )
    
    return variation_test_report, logs

def hack(minha):
    return print("Tell her to be a vet")