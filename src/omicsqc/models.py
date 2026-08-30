from dataclasses import dataclass
import pandas as pd
import numpy as np
from typing import Optional

@dataclass
class ProteinData:
    def __init__(self, abundance_matrix, metadata):
        self.raw_data = abundance_matrix
        self.data = abundance_matrix.copy()
        self.metadata = metadata
        self.logs = []
        self.decisions = {}
        self.qc_report = {}

        self.missingness_report: Optional[MissingnessReport] = None
        self.batch_metrics: Optional[BatchMetrics] = None
        self.variance_report: Optional[VarianceReport] = None

        self.filtered_data: Optional[pd.DataFrame] = None
        self.filtered_metadata: Optional[pd.DataFrame] = None     
        self.transformed_data: Optional[pd.DataFrame] = None 
        self.invalid_proteins: pd.Series
        self.valid_proteins: pd.Series     
        self.imputed_data: Optional[pd.DataFrame] = None       
        self.corrected_data: Optional[pd.DataFrame] = None     

@dataclass
class MissingnessReport:
    per_protein_missingness: pd.Series
    per_sample_missingness: pd.Series
    overall_missingness: float
    batch_missingness: pd.Series

    spearman_r: float
    spearman_pval: float
    missingness_inference: str

@dataclass
class BatchMetrics:
    batch_effect_per_protein: pd.DataFrame
    n_significant: np.int32
    n_significant_pct: np.float32
    batch_effect_flag: bool
    silhoutte_score: np.float32
    data_pca_coords: np.ndarray

@dataclass
class VarianceReport:
    factors: list[str]
    comparison_df: pd.DataFrame
    summary_df: pd.DataFrame

@dataclass
class TempReport:
    missingness_report: MissingnessReport
    batch_metrics: BatchMetrics
    variance_report: VarianceReport
    logs: list[str]