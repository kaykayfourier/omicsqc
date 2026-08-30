# omicsqc

## A Python package for automated quality control and preprocessing of mass spectrometry-based proteomics data.

Mass spectrometry proteomics studies measure thousands of proteins simultaneously across hundreds of samples. The raw data these experiments produce is systematically messy in ways specific to how the technology works — proteins below instrument detection thresholds are absent entirely, samples processed on different days carry technical artifacts that obscure real biology, and no standardized automated pipeline exists to handle these issues reproducibly from raw input to analysis-ready output.

omicsqc addresses this gap. It takes two standard input files — a protein intensity matrix and a sample metadata file — and produces a cleaned, normalized, batch-corrected matrix alongside a full audit report of every decision made during preprocessing. Every transformation is documented, every parameter is configurable, and every intermediate state is preserved so the user can inspect what changed at each stage.

The package is designed for proteomics researchers who want reproducible preprocessing without manually stitching together disconnected R and Python tools, and for computational biologists who need a clean, extensible foundation they can build domain-specific analyses on top of.

## What Has Been Completed

Data container — a ProteinData class that carries the expression matrix, metadata, all intermediate processing states, a decision log, and QC report objects as a single object throughout the pipeline.

Ingestion layer — loads and aligns the protein intensity matrix and metadata by sample ID, handles common ID column naming variations, replaces zeros with NaN, compresses memory via float32, validates for duplicate IDs, and returns a populated ProteinData object ready for QC.

QC layer — three independent, individually callable modules that together characterize data quality:

- Missingness analysis computes per-protein, per-sample, per-batch, and overall missingness rates, determines the statistical mechanism of missingness (MNAR, MAR, or MCAR) via Spearman correlation between detection rate and mean observed intensity, and outputs an imputation method recommendation grounded in that determination
- Batch metrics runs vectorized one-way ANOVA or Kruskal-Wallis across all proteins simultaneously to quantify what fraction of the proteome is significantly affected by batch, classifies batch effect severity as low, moderate, or severe, and stores per-protein F-statistics and significance flags
- Variance partitioning computes per-protein R² for any set of metadata factors (batch, diagnosis, sex, PMI) and produces a factor comparison showing which variables explain the most variance across the proteome — the key diagnostic for whether batch dominates over biology

All three run through a single compute_qc_pipeline() call or independently. Results attach directly to the ProteinData object and are accessible as typed dataclass attributes.

Preprocessing layer (partially complete) — filtering removes or flags proteins and samples exceeding configurable missingness thresholds, log2 transformation and median normalization are implemented and verified.

Remaining: imputation, batch correction, export, reporting, and Streamlit dashboard.