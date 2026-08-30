import numpy as np
import pandas as pd

from ..config import QCconfigs
from ..models import MissingnessReport

#Determining MAR vs MCAR vs MNAR nature of the dataset
#Returns spearman coefficient, spearman p-value, missingness mechanism and imputation_recommendation
def missingness_nature(data_df,  pval_threshold = QCconfigs.missingness_pval_threshold, corr_threshold = QCconfigs.missingness_corr_threshold, ):
    log2_df = np.log2(data_df.replace(0.0, np.nan))

    detection_rate = log2_df.notna().mean(axis = 1)
    mean_intensity_detected = log2_df.mean(axis=1, skipna=True)

    mask = detection_rate > 0
    detection_rate = detection_rate[mask]
    mean_intensity_detected = mean_intensity_detected[mask]
    from scipy.stats import spearmanr
    corr, pval = spearmanr(detection_rate, mean_intensity_detected)
    spearman_r = corr
    spearman_pval = pval
     
    if corr > corr_threshold and pval < pval_threshold:
        inference = "Evidence of intensity-dependent missingness, consistent with left-censoring  "
    elif abs(corr) < 0.1:
        inference = "Extremely poor relation between protein intensity and detection rate, consistent future reportings may infer MCAR"
    else:
        inference = "Unjustifiable evidence of intensity-dependent missingness, missiningness may be of MAR nature"

    return spearman_r, spearman_pval, inference

# Complete Function that computes all variables of the Missingness Report and returns dataclass of structure specified in models.py

def compute_missingness(data_df, meta_df, pval_threshold = QCconfigs.missingness_pval_threshold, corr_threshold = QCconfigs.missingness_corr_threshold, ):
    
    overall_missingness = data_df.isna().mean().mean() * 100
    #print(overall_missingness)
    if overall_missingness < 1e-8:    
        logs = "missingness analysis not applicable because no missing values remain"
        report = None
        return report, logs
    
    else:
        per_protein_missingness = pd.Series(data_df.isna().mean(axis=1) * 100)
        per_sample_missingness = pd.Series(data_df.isna().mean(axis = 0) * 100)
        possible_batch_cols = ["Batch_overall","batch","Batch","Batches","MS_Batch","RunBatch"]
    
        batch_col = next(
            (c for c in meta_df.columns
            if c in possible_batch_cols or "batch" in c.lower()),
            None
        )

        sample_to_batch = meta_df[batch_col].dropna()
        common = data_df.columns.intersection(sample_to_batch.index)
        expr_aligned = data_df[common]

        batch_missingness = (expr_aligned.isna()
                            .T
                            .groupby(sample_to_batch[common].values)
                            .mean()
                            .mean(axis=1) * 100).sort_values(ascending = False)
            
            
        batch_missingness = pd.Series(batch_missingness).sort_values(ascending = False)
        spearman_r, spearman_pval, inference = missingness_nature(data_df, pval_threshold = pval_threshold, corr_threshold = corr_threshold)
        
        report = MissingnessReport(per_protein_missingness, per_sample_missingness, overall_missingness, batch_missingness, spearman_r, spearman_pval, missingness_inference= inference)
        logs = (
            f"""
        [Missingness Analysis]
        Per-protein missingness (mean): {per_protein_missingness.mean():.2f}%
        Per-sample missingness (mean):  {per_sample_missingness.mean():.2f}%
        Overall missingness:            {overall_missingness:.2f}%
        Spearman coefficient:           {spearman_r:.4f}
        Spearman p-value:               {spearman_pval:.4e}
        Inferred missingness nature: {inference}
        """.strip()
        )
        
        return report, logs
