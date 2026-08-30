import numpy as np
import pandas as pd
from ..ingestion import ProteinData
from omicsqc.qc_layer.missingness_analysis import compute_missingness
from omicsqc.qc_layer.batch_metrics import batch_variation_test
from omicsqc.qc_layer.variation_test import compare_explained_variances
from ..config import QCconfigs
from ..models import TempReport


def compute_qc_pipeline(protdata: ProteinData, data_df:pd.DataFrame, meta_df : pd.DataFrame, protein_threshold_pct: np.float32, variation_factors: list[str],
                         qc_configs: QCconfigs = None) -> ProteinData:
    
    if qc_configs == None:
        qc_configs = QCconfigs

    logs = []
    protdata.missingness_report, log = compute_missingness(data_df= data_df,
                                            meta_df= meta_df,
                                            pval_threshold=qc_configs.missingness_pval_threshold,
                                            corr_threshold= qc_configs.missingness_corr_threshold,  )
    logs.append(log)
    
    protdata.batch_metrics, log = batch_variation_test(data_df= data_df, meta_df=meta_df, protein_threshold_pct= protein_threshold_pct,
                                         batch_effect_method= qc_configs.batch_effect_method,
                                         alpha= qc_configs.alpha)
    logs.append(log)
    
    protdata.variation_test_report, log = compare_explained_variances(data_df= data_df, meta_df=meta_df, factors= variation_factors)
    logs.append(log)
    protdata.logs = logs
    return protdata


def compute_qc_individual(data_df:pd.DataFrame, meta_df : pd.DataFrame, protein_threshold_pct: np.float32, variation_factors: list[str],
                         qc_configs: QCconfigs = None):
        
    if qc_configs == None:
        qc_configs = QCconfigs

    logs = []

    missingness_report, log = compute_missingness(data_df, meta_df,
                                            pval_threshold=qc_configs.missingness_pval_threshold,
                                            corr_threshold= qc_configs.missingness_corr_threshold,  )
    logs.append(log)

    batch_metrics, log = batch_variation_test(data_df, meta_df, protein_threshold_pct= protein_threshold_pct,
                                         batch_effect_method= qc_configs.batch_effect_method,
                                         alpha= qc_configs.alpha)
    logs.append(log)

    variation_test_report, log = compare_explained_variances(data_df, meta_df, factors= variation_factors)
    logs.append(log)

    temp_report = TempReport(missingness_report= missingness_report, batch_metrics= batch_metrics, variance_report= variation_test_report, logs= logs)

    return temp_report

