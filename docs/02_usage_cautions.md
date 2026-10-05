# omicsqc: Usage Notes, Assumptions and Known Limitations

This document lists what a user needs to know before running omicsqc and trusting its output. It covers the assumptions the code makes, the settings that are judgement calls, and the known weaknesses.

## 1. Input data assumptions

- **Only two inputs:** an expression matrix (proteins as rows, samples as columns) and a metadata table (one row per sample).
- **Only the TMT ROSMAP/Banner format is supported.** Sample IDs must look like `rmap_b32.129N` and the batch column must hold values like `rmap_b32`. Anything else fails at ingestion with a format error. This is on purpose: guessing other naming schemes is how wrong batch labels slip in.
- **You must give the real technical batch column.** In the original metadata, the column called `batch` only says `rmap` or `bnnr`, which is the cohort, not the batch. Using it would treat 67 plexes as 2 groups. Create a column holding the part of the sample ID before the dot and pass that as `batch_col`.
- **The expression matrix must be raw, unnormalised intensities.** The transform step refuses data that is already log2-transformed. Zeros in the expression matrix are treated as missing values. Zeros in the metadata are left alone, since 0 is a real value there (for example sex coded 0/1).
- **If the expression matrix has no protein ID index, the first column is used as the protein ID.**

## 2. Which data goes into which function

Getting this wrong is the easiest way to get misleading numbers.

| Function | Expects |
|---|---|
| `log2_median_centering` | Raw intensities |
| `compute_missingness` | Raw intensities (it applies log2 itself) |
| `batch_variation_test`, `compare_explained_variances` | log2-transformed and median-centred data |
| `filter_data` | log2-transformed and median-centred data |
| `impute_*` | Filtered, log2-transformed, median-centred data |

The statistical functions raise an error if the maximum value is above 100 (that can't be log2 data), and warn if sample medians aren't close to 0.

**Run QC on the filtered data with its NaNs, not on the imputed data.** MinProb fills missing plexes with random low values that are much noisier than real values within a plex. On this dataset that lowered the median batch eta² of the affected proteins from 0.976 to 0.914, which makes the batch effect look smaller than it is. The imputed matrix is only an input for batch correction methods that can't handle missing values.

## 3. Settings that are judgement calls

These defaults were chosen by looking at the real data, but they are not universal rules.

| Setting | Default | What it means | When to change it |
|---|---|---|---|
| `max_protein_missing_pct` (filtering) | 80% | On this data: protein must be seen in about 14 of 67 plexes | Use 50% if you want less imputed data (keeps 9,304 instead of 10,448 proteins) |
| `max_sample_missing_pct` (filtering) | 50% | Removes very poor samples | No sample was removed on this data |
| `min_batches` (ANOVA) | 3 | Proteins seen in fewer plexes aren't tested | Proteins seen in only 3 plexes can still get extreme F values. Raise it (for example to 25% of the number of batches) for stricter results. |
| `eta_threshold` (batch flag) | 0.1 | A protein is flagged if batch explains at least 10% of its variance | Before correction almost every protein is above this. It becomes useful when comparing before and after correction. |
| `min_obs` (variance explained) | 30 | Proteins with fewer observed samples get no R² | Lower values let sparse proteins produce unstable, extreme values |
| `detection_cutoff` (mixed imputation) | 0.9 | Proteins seen in under 90% of samples are treated as "too low to detect" | Based on the detection curve, which reaches ~95-99% detection at higher intensities |
| `corr_threshold` (missingness) | 0.3 | Spearman rho above this is called intensity-dependent | |

## 4. How to read the results

**Batch test**
- With 558 samples, even tiny effects give tiny p-values. Look at **eta²** (how big the effect is), not just whether it is significant.
- A p-value of exactly 0 means it is too small for the computer to store (below about 1e-308). Rank those proteins by F or eta² instead.
- Proteins marked `low_coverage` were not tested, and their results are NaN on purpose.

**Batch vs cohort**
- Every plex belongs to only one cohort (ROSMAP or Banner). So a large batch effect can partly be a cohort difference (different biology, sample handling or processing), not purely the instrument run.
- To separate them, run the test inside each cohort separately and compare. Signal that stays inside a cohort is batch.

**Missingness**
- The test can only say whether missingness **depends on intensity** (yes on this data, rho = 0.88). It **cannot** say whether data is MAR or MNAR, because those two can't be told apart from the observed data alone. The verdict text says this.
- On this data missingness is plex-level: a protein is either measured in every sample of a plex or in none. Because of that:
  - The batch row in the covariate table (effect size 0.9998) is true by construction and tells you nothing new.
  - Sample-level p-values for other covariates are inflated, because the 8-9 samples in a plex share one outcome. The real number of independent units is 67 plexes, not 558 samples. Read `effect_size`, not the p-value.
- The mean intensity of a rarely-detected protein is biased upward, because only its brighter measurements were seen. This makes the Spearman test somewhat conservative.

**Variance explained**
- The value is **adjusted R²**. 0 means "no better than chance", and small negative values are normal.
- Values for different factors **must not be added up**: the factors overlap (batch is inside cohort, and diagnosis may be unevenly spread over plexes).
- `cogdx` is only recorded for ROSMAP samples. With `complete_cases=True` (default), including it restricts *every* factor to the 360 ROSMAP samples (45 plexes), and the batch number then describes ROSMAP only.
- Before correction, a biological factor's R² can include some batch signal. After correction, expect batch R² to drop and biological R² to **rise**, since the same biological effect becomes a bigger share of the remaining variance. If biological R² falls, the correction may have removed real signal.

**Imputation**
- MinProb values are random. Pass a `seed` to get the same result every run.
- Imputed values are guesses. For proteins in the `left_censored` group, up to 80% of values can be imputed.
- The kNN step runs on thousands of proteins and can take about a minute.

## 5. Statistical assumptions

- **ANOVA** assumes the data in each batch is roughly normal with similar spread. log2 + centering helps a lot, but with 8-9 samples per batch this can't really be checked per protein. Kruskal-Wallis is available if this is a concern.
- **Detection of "already transformed" data is a rule of thumb** (negative values or max < 50). It works for this dataset's raw scale (values up to ~1e10) but is not a guarantee for every dataset.
- **Imputation choices assume** that rarely-detected proteins are missing because they were too low (left-censored), and that occasional gaps are random. These are standard assumptions in proteomics but cannot be proven from the data.

## 6. Known limitations and unfinished parts

- Only the TMT ROSMAP/Banner naming convention is supported.
- Batch correction is not implemented yet, so there is no "after correction" result.
- The pipeline wrapper (`qc_metrics.py`) is outdated and will fail if called. Use the module functions directly for now.
- `config.py` has two settings written with trailing commas (`MinProb_lmbda`, `MinProb_delta`), which makes them tuples. They are not used by the new imputation code, but should be removed. `missingness_pval_threshold` is also unused.
- Kruskal-Wallis has no effect size and ignores `min_batches`.
- `min_batches` is not exposed in `batch_variation_test`; to change it you currently have to call `vectorized_anova` directly.
- The ingestion order checks batch sizes before checking the naming format, so a badly formatted batch column can show a batch-count error instead of the clearer format error.
- In a Jupyter notebook, changes to package files are not picked up until the kernel restarts (or `%autoreload 2` is used). This caused confusing "unexpected keyword argument" errors during development.
