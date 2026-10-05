# omicsqc: Architecture and Modules

*Status as of 5 October 2026. Batch correction is not implemented yet.*

## 1. What the package does

omicsqc is a Python package for quality control (QC) and preprocessing of mass-spectrometry proteomics data, mainly to find and later remove **batch effects**.

A batch effect is a technical difference between groups of samples that were processed together. In TMT proteomics, each **plex** (a group of 8-9 samples measured in one instrument run) is a batch. Samples in the same plex end up looking more alike than samples from different plexes, even when the biology is the same. If this is not handled, the batch can hide or fake real biological differences.

The package takes two files and does the following:

1. Loads and checks the data (ingestion).
2. Measures how bad the batch effect and missing data are (QC layer).
3. Cleans the data: log2 transform, filtering and imputation (preprocessing layer).
4. Corrects the batch effect (batch layer, not built yet).
5. Runs the QC layer again so the before and after numbers can be compared.

The main idea is that these steps, which are usually spread across many notebooks and tools, are organised into layers with checks and a record of what was done at each step.

## 2. Input data

The package currently supports one data format: the **TMT ROSMAP and Banner** dataset ([Synapse syn52854669](https://www.synapse.org/Synapse:syn52854669)).

| File | Shape | Notes |
|---|---|---|
| Expression matrix (`ProteinAbundances.csv`) | proteins x samples | Raw, unnormalised FragPipe abundances. 11,958 proteins x 558 samples. |
| Metadata (`meta_data.csv`) | samples x traits | One row per sample: batch, diagnosis, sex, APOE genotype and so on. |

Sample IDs follow `<cohort>_<batch>.<channel>`, for example `rmap_b32.129N`:
- `rmap` / `bnnr` is the cohort (ROSMAP or Banner),
- `b32` is the plex number,
- `129N` is the TMT channel inside that plex.

So the real technical batch of a sample is the part before the dot (`rmap_b32`). The dataset has 67 plexes (45 ROSMAP, 22 Banner) with 8-9 samples each.

## 3. Package layout

```
src/omicsqc/
├── ingestion.py              loading, alignment and format checks
├── models.py                 ProteinData container and report dataclasses
├── config.py                 default parameters (QCconfigs)
├── qc_layer/
│   ├── missingness_analysis.py   how much data is missing and in what pattern
│   ├── batch_metrics.py          per-protein ANOVA / Kruskal-Wallis across batches
│   ├── variation_test.py         variance explained (adjusted R²) per factor
│   └── qc_metrics.py             old pipeline wrapper (outdated, see section 8)
├── preprocessing/
│   ├── transform.py              log2 + median centering
│   ├── filtering.py              drop samples/proteins with too much missing data
│   └── impute.py                 MinProb, kNN and mixed imputation
└── batch/
    └── batch_correction.py       empty for now
tests/                            pytest files for every module
```

## 4. Workflow

```
raw data ──ingest──► ProteinData
                         │
                    log2 + median centering ──► transformed_data
                         │
                    QC run 1: "diagnostic"   (on transformed data)
                         │
                    filtering ──► filtered_data
                         │
                    QC run 2: "pre-correction"  (on filtered data, with NaNs)
                         │
                    imputation ──► imputed_data   (only because correction needs it)
                         │
                    batch correction ──► corrected_data   (not built yet)
                         │
                    QC run 3: "corrected"
```

The comparison that shows whether correction worked is **QC run 2 vs QC run 3**. Both are computed on the same proteins and samples, so the numbers are directly comparable. QC run 1 is only there to understand the data and choose the filter and imputation settings.

## 5. Two ways to use it

**Fragmented mode (the core functions).** Every module function takes plain pandas DataFrames or Series and returns a result plus a log string. It does not need `ProteinData`. This is for experimenting, for example trying several filter thresholds.

**Pipeline mode (with `ProteinData`).** `ProteinData` stores each stage of the data and the QC reports, so the whole run is traceable. The pipeline wrappers only read a stage, call the same core function and store the result, so both modes always give the same numbers. These wrappers are still to be written (section 8).

## 6. The data container: `ProteinData` (models.py)

`ProteinData` holds the data and everything computed from it.

| Attribute | Meaning |
|---|---|
| `raw_data` (also `data`) | Expression matrix after ingestion |
| `transformed_data`, `filtered_data`, `imputed_data`, `corrected_data` | One matrix per stage. `None` means "not computed yet". |
| `metadata`, `batch_col` | Metadata and the name of the technical batch column |
| `valid_proteins`, `invalid_proteins` | Result of filtering |
| `missingness_report_<run>`, `batch_metrics_<run>`, `variance_report_<run>` | QC reports for each run (`diagnostic`, `imputed`, `batch_corrected`) |
| `logs`, `decisions` | Text logs and choices made |

| Method | What it does |
|---|---|
| `get_stage(stage)` | Returns `(matrix, metadata)` for a stage. The metadata is cut to the same samples, so the two can never be misaligned. Raises an error if the stage hasn't been computed. |
| `set_stage(stage, df, overwrite=False)` | The only way to save a stage. Refuses to overwrite an existing stage (this stops, for example, log2-transforming twice) and checks all samples exist in the metadata. |
| `set_qc(name, missingness, batch_metrics, variance)` | Saves the reports of one QC run. |

Report dataclasses: `MissingnessReport`, `BatchMetrics` and `VarianceReport` (fields listed with each module below).

## 7. Modules

### 7.1 Ingestion: `ingest_data` (ingestion.py)

```python
ingest_data(protein_data, meta_df, batch_col, sample_id_indexed=False,
            convention="TMT", float32_compress=True) -> ProteinData
```

| Parameter | Meaning |
|---|---|
| `protein_data` | Expression matrix (proteins x samples). If it has no protein index, the first column is used. |
| `meta_df` | Metadata table |
| `batch_col` | Name of the metadata column holding the technical batch (for example `rmap_b32`) |
| `sample_id_indexed` | `True` if the metadata index already holds sample IDs. If `False`, the column whose values best match the expression columns is used. |
| `convention` | Naming convention to validate. Only `"TMT"` exists. |
| `float32_compress` | Store the matrix as float32 to halve memory use |

Steps: set protein index → find the sample ID column → check `batch_col` exists → check for duplicate IDs → keep only samples present in both files (in the same order) → check there are at least 2 batches with at least 2 samples each → validate the TMT naming → replace zeros with NaN in the expression matrix only → compress.

TMT validation (`TMT_convention_validator`) checks that every batch value looks like `<cohort>_b<number>`, that every sample ID looks like `<cohort>_b<number>.<channel>`, and that the batch part of each sample ID equals that sample's batch value.

**Output:** a `ProteinData` object.

### 7.2 Transform: `log2_median_centering` (preprocessing/transform.py)

```python
log2_median_centering(df) -> DataFrame
```

Takes raw intensities, applies log2, then subtracts each sample's median so all samples are centred at 0. NaNs stay NaN, and zeros become NaN with a warning.

Guards: it raises an error if the data has negative values or a maximum below 50, because that means it was already transformed. This prevents transforming twice.

Why: the statistical tests in the QC layer assume roughly normal data with similar spread. Raw intensities are heavily skewed, and log2 fixes most of that. Median centering removes differences in total sample loading.

### 7.3 Missingness analysis: `compute_missingness` (qc_layer/missingness_analysis.py)

```python
compute_missingness(data_df, meta_df, batch_col, covariates=None,
                    corr_threshold=0.3) -> (MissingnessReport or None, logs)
```

| Parameter | Meaning |
|---|---|
| `data_df` | **Raw** intensities (the function applies log2 itself) |
| `batch_col` | Batch column name |
| `covariates` | Metadata columns to test against missingness. Default: just the batch column. |
| `corr_threshold` | Spearman rho above which missingness is called intensity-dependent |

What it computes:
- Missing % per protein, per sample, per batch and overall.
- **Intensity dependence:** Spearman correlation between how often a protein is detected (per plex) and its mean intensity, plus a **detection curve** (detection rate in each intensity decile). If low-intensity proteins are missing more often, the data is "left-censored": values are missing because they were below the instrument's detection limit.
- **Partial fraction:** the share of protein-plex cells that are only partly missing. Near 0 means proteins go missing for whole plexes at once.
- **Covariate table:** whether sample missingness differs by a metadata column (Kruskal-Wallis for categories, Spearman for numbers).

`MissingnessReport` fields: `per_protein_missingness`, `per_sample_missingness`, `overall_missingness`, `batch_missingness` (all in %), `spearman_r`, `spearman_pval`, `missingness_inference` (text verdict), `detection_curve` (0-1), `n_proteins_tested`, `partial_fraction` (0-1), `covariate_table`.

Returns `None` if nothing is missing.

### 7.4 Batch test: `batch_variation_test` (qc_layer/batch_metrics.py)

```python
batch_variation_test(data_df, meta_df, batch_col, batch_effect_method="anova",
                     alpha=0.05, eta_threshold=0.1) -> (BatchMetrics, logs)
```

Input must be log2-transformed and median-centred.

For each protein it asks: **does the mean level differ between batches?**

- `"anova"` uses `vectorized_anova`: a one-way ANOVA computed for all proteins at once with NumPy. For each protein it splits the variation into *between batches* and *within batches* and computes F = (between variation per degree of freedom) / (within variation per degree of freedom). It also gives eta² = share of the protein's variance explained by batch (0 to 1).
- `"kruskal_wallis"` uses a rank-based test that doesn't assume normality (slower, no eta²).

Details:
- Degrees of freedom are counted per protein, using only the batches where that protein was observed.
- Proteins seen in fewer than `min_batches=3` batches are not tested (`low_coverage=True`).
- p-values are corrected for testing thousands of proteins using **Benjamini-Hochberg** (q-values).
- `significant` = q < alpha. `Flag` = significant **and** eta² ≥ `eta_threshold`, which means the effect is both real and big.

Per-protein table columns: `F`, `p_value`, `q_value`, `eta_sq`, `n_batches`, `low_coverage`, `significant`, `Flag`.

`BatchMetrics` fields: `batch_effect_per_protein` (the table), `n_testable`, `n_significant`, `n_significant_pct`, `n_flagged`, `n_flagged_pct`. Percentages are out of testable proteins.

### 7.5 Variance explained: `compare_explained_variances` (qc_layer/variation_test.py)

```python
compare_explained_variances(data_df, meta_df, factors, min_groups=2,
                            min_obs=30, complete_cases=True) -> (VarianceReport, logs)
```

For each factor (batch, diagnosis, sex and so on) and each protein, it computes the **adjusted R²**: how much of the protein's variance that factor explains, corrected for the number of groups. This lets batch be compared with biology on the same scale.

| Parameter | Meaning |
|---|---|
| `factors` | Categorical metadata columns to compare |
| `min_groups` | Minimum groups a protein must be seen in |
| `min_obs` | Minimum observed samples per protein |
| `complete_cases` | Use only samples that have a value for every factor, so all factors are measured on the same samples |

`VarianceReport` fields: `factors`, `comparison_df` (proteins x factors), `summary_df` (per factor: `n_groups`, `n_samples`, `n_proteins_tested`, mean/median/std/quantiles of adjusted R²).

### 7.6 Filtering: `filter_data` (preprocessing/filtering.py)

```python
filter_data(df, max_protein_missing_pct=80.0, max_sample_missing_pct=50.0,
            batch=None) -> (filtered_df, protein_mask, sample_mask, logs)
```

Removes samples with more than 50% missing values first, then proteins with more than 80% missing (computed on the remaining samples). Proteins with no values at all are always removed. If `batch` is given, it warns when a batch loses all its samples or is left with fewer than 2.

The masks are True/False per protein and per sample (True = kept), so it's always clear what was removed.

### 7.7 Imputation (preprocessing/impute.py)

All three take the filtered, log2, centred matrix and return a matrix with no NaNs.

| Function | Use | How it works |
|---|---|---|
| `impute_minprob(df, q=0.01, tune_sigma=1.0, seed=None)` | Values missing because they were too low to detect | Fills each sample's gaps with random values drawn near that sample's 1st percentile |
| `impute_knn(df, k=10)` | Values missing at random | Fills a gap using the average of the 10 most similar proteins in the same sample |
| `impute_mixed(df, detection_cutoff=0.9, ...)` | Default | Proteins seen in less than 90% of samples → MinProb; the rest with gaps → kNN; complete proteins untouched. Returns `(imputed, groups, logs)`, where `groups` labels each protein `complete`, `left_censored` or `random`. |

## 8. Current status

| Part | Status |
|---|---|
| Ingestion | Done and tested |
| QC layer (missingness, batch test, variance) | Done and tested |
| Preprocessing (transform, filter, impute) | Done and tested |
| `ProteinData` stages | Done |
| Pipeline wrappers | Not yet. `qc_metrics.py` still calls old function signatures and would fail if run. To be rewritten using `get_stage` / `set_stage`. |
| Batch correction | Not started. Options under consideration: per-plex median centering (handles missing values) and ComBat via InMoose (needs complete data). |
| Reporting / dashboard | Not started |

## 9. Results on the real dataset (diagnostic QC)

| Measure | Result |
|---|---|
| Samples / proteins / plexes | 558 / 11,958 / 67 |
| Overall missing values | 22.3% |
| Detection vs intensity (Spearman rho) | 0.88: low-abundance proteins go missing much more often |
| Partly-missing protein-plex cells | 0%: proteins are missing for whole plexes at a time |
| Batch effect, median variance explained by batch (adjusted R², ROSMAP samples) | 0.95 |
| Variance explained by diagnosis, sex, APOE (median) | about 0 |
| After filtering (80% default) | 10,448 proteins kept, missing values down to 12.3% |

So batch explains about 95% of a typical protein's variance before correction, which is the main reason batch correction is needed.
