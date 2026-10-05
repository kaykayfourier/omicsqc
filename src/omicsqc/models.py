from dataclasses import dataclass
import pandas as pd
import numpy as np
from typing import Optional

STAGES = {
    "raw":         "raw_data",
    "transformed": "transformed_data",
    "filtered":    "filtered_data",
    "imputed":     "imputed_data",
    "corrected":   "corrected_data",
}

QC_STAGES = {
    "diagnostic":      "transformed",
    "imputed":         "imputed",
    "batch_corrected": "corrected",
}



class ProteinData:
    def __init__(self, abundance_matrix: pd.DataFrame, metadata: pd.DataFrame, batch_col: str):
        self.metadata = metadata
        self.batch_col = batch_col
        self.logs = []
        self.decisions = {}
        self.qc_report = {}
 
        # one matrix per pipeline stage. None means "not computed yet".
        self.raw_data: pd.DataFrame = abundance_matrix
        self.transformed_data: Optional[pd.DataFrame] = None
        self.filtered_data: Optional[pd.DataFrame] = None
        self.imputed_data: Optional[pd.DataFrame] = None
        self.corrected_data: Optional[pd.DataFrame] = None
 
        # filter results
        self.valid_proteins: Optional[pd.Series] = None
        self.invalid_proteins: Optional[pd.Series] = None
 
        # QC reports, one set per QC run (see QC_STAGES)
        self.missingness_report_diagnostic: Optional[MissingnessReport] = None
        self.batch_metrics_diagnostic: Optional[BatchMetrics] = None
        self.variance_report_diagnostic: Optional[VarianceReport] = None
 
        self.missingness_report_imputed: Optional[MissingnessReport] = None
        self.batch_metrics_imputed: Optional[BatchMetrics] = None
        self.variance_report_imputed: Optional[VarianceReport] = None
 
        self.missingness_report_batch_corrected: Optional[MissingnessReport] = None
        self.batch_metrics_batch_corrected: Optional[BatchMetrics] = None
        self.variance_report_batch_corrected: Optional[VarianceReport] = None
 
    @property
    def data(self) -> pd.DataFrame:
        return self.raw_data
 
    @property
    def filtered_metadata(self) -> Optional[pd.DataFrame]:
        if self.filtered_data is None:
            return None
        return self.metadata.loc[self.filtered_data.columns]
 
    def get_stage(self, stage: str):
        """Returns (matrix, metadata) for a stage. The metadata is subset to the matrix's samples."""
        if stage not in STAGES:
            raise ValueError(f"Unknown stage '{stage}'. Use one of {list(STAGES)}.")
        df = getattr(self, STAGES[stage])
        if df is None:
            raise ValueError(f"Stage '{stage}' hasn't been computed yet.")
        return df, self.metadata.loc[df.columns]
 
    def set_stage(self, stage: str, df: pd.DataFrame, overwrite: bool = False):
        """The only way a pipeline step should write a stage."""
        if stage not in STAGES:
            raise ValueError(f"Unknown stage '{stage}'. Use one of {list(STAGES)}.")
        attr = STAGES[stage]
        if getattr(self, attr) is not None and not overwrite:
            raise ValueError(f"Stage '{stage}' already exists. Pass overwrite=True to replace it.")
        unknown = df.columns.difference(self.metadata.index)
        if len(unknown):
            raise ValueError(f"{len(unknown)} samples in the '{stage}' matrix are not in the metadata.")
        setattr(self, attr, df)
 
    def set_qc(self, name: str, missingness=None, batch_metrics=None, variance=None):
        """Stores the QC reports of one run under its name (diagnostic, imputed, batch_corrected)."""
        if name not in QC_STAGES:
            raise ValueError(f"Unknown QC run '{name}'. Use one of {list(QC_STAGES)}.")
        setattr(self, f"missingness_report_{name}", missingness)
        setattr(self, f"batch_metrics_{name}", batch_metrics)
        setattr(self, f"variance_report_{name}", variance)
    

@dataclass
class MissingnessReport:
    per_protein_missingness: pd.Series
    per_sample_missingness: pd.Series
    overall_missingness: float
    batch_missingness: pd.Series
    spearman_r: float
    spearman_pval: float
    missingness_inference: str
    detection_curve: pd.Series
    n_proteins_tested: int
    partial_fraction: float
    covariate_table: pd.DataFrame

@dataclass
class BatchMetrics:
    batch_effect_per_protein: pd.DataFrame
    n_testable: int
    n_significant : int
    n_significant_pct: float
    n_flagged: Optional[int]
    n_flagged_pct: float   # NaN when n_flagged is None
    
    
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