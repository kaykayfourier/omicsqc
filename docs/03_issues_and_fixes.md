# omicsqc: Problems Encountered and How They Were Solved

This is a log of the main problems I ran into while building omicsqc, what caused them, and what I changed. They are grouped by part of the package, roughly in the order I worked on them.

## 1. Ingestion

**1.1 The "batch" column was not the real batch.**
The metadata column `batch` only contains `rmap` or `bnnr`, which are the two cohorts. The real batch (plex) is only visible inside the sample ID: in `rmap_b32.129N`, the batch is `rmap_b32`. All my batch tests were therefore comparing 2 cohorts instead of 67 plexes.
*Fix:* The user must now pass an explicit technical batch column. Ingestion checks that every batch value has the form `<cohort>_b<number>` and that it matches the batch part of the sample ID. I decided not to auto-detect batches from sample IDs, since every dataset names samples differently and a wrong guess fails silently.

**1.2 Naming conventions differ between datasets.**
I couldn't find a standard naming scheme for proteomics sample IDs. Only the TMT channel part (`127N`, `129C` and so on) is standard.
*Fix:* I limited the package to the TMT ROSMAP/Banner format and wrote that in the assumptions. Conventions are selected by name (`convention="TMT"`), so another one can be added later as one extra validator.

**1.3 Finding the sample ID column relied on guessing names.**
The old code looked for column names containing "id" or "sample", which could pick the wrong column (for example `individualID`).
*Fix:* Each candidate column is scored by how many of its values match the expression matrix's column names, and the best match is used.

**1.4 Several small bugs** found by reviewing the code and writing tests:
- `max(scores, key=scores.get())` called the function instead of passing it.
- The duplicate-sample check looked at protein IDs (rows) instead of sample IDs (columns).
- float32 compression was applied to a variable that was never returned.
- Common samples were taken from a Python `set`, which has no fixed order, so data and metadata could end up in different orders. Now the expression matrix's order is kept and metadata is reordered to match.
- Zeros were replaced by NaN in the **metadata** too, which deleted real values such as sex coded 0/1. Now this only happens in the expression matrix.

**1.5 Default setting caused "No common Sample IDs" on real data.**
`sample_id_indexed=True` by default meant the metadata's 0..N row numbers were treated as sample IDs.
*Fix:* Default changed to `False`, so the ID column is detected automatically.

I wrote 14 pytest cases that each corrupt a small toy dataset in one way (duplicate IDs, wrong format, single batch and so on) and check the right error appears. All pass, and the real data ingests 558 of 558 samples.

## 2. Batch test (ANOVA)

**2.1 ANOVA was run on raw intensities.**
ANOVA assumes roughly normal data with similar spread in each group. Raw intensities are very skewed.
*Fix:* Data is log2-transformed and median-centred first. The ANOVA now raises an error if the data looks raw (max > 100) and warns if samples aren't centred.

**2.2 Wrong degrees of freedom.**
The code used the total number of batches (k) for every protein. But a protein missing from some batches only has data in fewer groups, so its degrees of freedom must be smaller. Using the global k gave wrong p-values, especially for proteins with batch-dependent missing data.
*Fix:* k is now counted per protein (only batches where it was observed).

**2.3 The "Flag" column had no statistical meaning.**
It flagged proteins with F above the average F and an uncorrected p < 0.05, while the `significant` column next to it used Bonferroni. Two different definitions of "affected" in one table.
*Fix:* Bonferroni was replaced with **Benjamini-Hochberg** (it controls the share of false hits instead of the chance of any false hit, which suits screening thousands of proteins). I added **eta²** as an effect size. Now `Flag` = significant and eta² ≥ threshold.

**2.4 A bug that only showed up with more than one protein.**
`np.nansum(group)` was missing `axis=1`, so it summed the whole batch across all proteins instead of per protein. My first toy test passed by luck, because of how the numbers lined up.
*How it was found:* a second toy protein where every batch mean is 0 should have given F = 0 but gave F = 1.45.
*Fix:* Added `axis=1`, and wrote 17 tests including a comparison against SciPy's `f_oneway` on data with missing values, plus tests that a protein's result doesn't change when other proteins are added or rows are reordered.

**2.5 The batch effect on real data looked "too strong".**
Top proteins had eta² ≈ 0.996, meaning almost all variance was between batches. I suspected a bug.
*What I found:* Looking at one protein (ZBTB2), values within a plex varied by ~0.2 log2 units, while the two plexes differed by ~5 log2 units (~35×). TMT channels in a plex are measured in the same run, so within-plex noise is small and plex-to-plex differences are large. The result was real, not a bug.

**2.6 Proteins seen in only 2-3 plexes topped the ranking.**
With very few samples, F can be extremely large and doesn't describe the 67-plex design.
*Fix:* Added `min_batches` (default 3). Proteins below it aren't tested and are marked `low_coverage`; how many were excluded is logged.

## 3. Missingness analysis

**3.1 The old function claimed to identify MAR / MNAR / MCAR.**
It correlated detection rate with mean intensity and labelled the result as one of the three types. But MAR and MNAR can't be told apart from observed data alone, so the "MAR" label could never be supported. Also, for proteins that are rarely detected only their brightest values are seen, which biases their mean upwards and weakens the correlation.
*Fix:* The function now only reports what it can show: whether missingness depends on intensity (and therefore rules out completely random missingness), plus a detection curve by intensity decile. The verdict text says MAR vs MNAR can't be distinguished. I looked at detection-probability modelling, but it is a research problem of its own and out of scope.

**3.2 Detection was counted per sample instead of per plex.**
In TMT, a protein is usually either identified in a whole plex or not at all. So detection is now counted per plex. On the real data, 0% of protein-plex cells are partly missing, which confirmed this.

**3.3 Testing missingness against batch gives a meaningless result.**
Since all samples in a plex share the same missingness, batch "explains" it perfectly (effect size 0.9998), and p-values are inflated because 8-9 samples count as 8-9 results when they're really one.
*Decision:* I tried a version that tests one value per plex, but it made the code much larger. I kept the simpler version and documented the limitation instead.

## 4. Variance explained

**4.1 R² was biased towards factors with many groups.**
Plain R² goes up with the number of groups even for pure noise (about (k−1)/(n−1)), so 67-level batch would beat 2-level diagnosis even if both had no real effect.
*Fix:* Switched to **adjusted R²**. Also reused the per-protein group count from the ANOVA, dropped samples with missing labels before computing the overall mean, and rejected numeric factors such as age (which would have created one group per value).

**4.2 `cogdx` showed impossible values (R² up to 0.93).**
*Cause 1:* `cogdx` is missing for all 198 Banner samples, so it was measured on different samples than the other factors.
*Cause 2:* Some proteins were observed in only ~25 samples across 3 plexes, so plex offsets leaked into any factor that happened to line up with them.
*Fix:* `complete_cases=True` (all factors use the same samples) and `min_obs=30`. After this, the median R² for diagnosis, sex and APOE was about 0 and batch was 0.95.

## 5. Preprocessing

**5.1 Filtering was far too strict.**
The old filter dropped any protein that was fully missing in at least one batch. Because missingness is plex-level, this removed almost every protein with *any* missing value, leaving nothing to impute.
*Fix:* That rule was removed. Samples are filtered first, then proteins by overall missingness. I tested thresholds from 20% to 95% on the real data. There was no clear "best" point, so I chose 80% as a lenient default (keeps 10,448 of 11,958 proteins, missing values 22.3% → 12.3%, all 67 plexes intact).

**5.2 Transform and filter were mixed into one function.**
`transform_data` did both, even though they happen at different points in the workflow, and could be run twice on the same data.
*Fix:* Split into `log2_median_centering` and `filter_data`. The transform now refuses data that already looks transformed. `ProteinData.set_stage` also refuses to overwrite a stage, which is a second guard.

**5.3 The old per-batch MinProb never actually used its per-batch logic.**
Since a missing protein is missing in its entire plex, the per-batch mean was always NaN and the code always fell back to global values.
*Fix:* Replaced with the standard per-sample MinProb, added kNN (neighbours are similar proteins), and a mixed function that sends rarely-detected proteins to MinProb and occasional gaps to kNN.

**5.4 Imputation changed the batch numbers.**
I compared the batch test before and after imputation. kNN didn't change anything, but for MinProb proteins eta² fell from 0.976 to 0.914, because random fills add noise inside each plex.
*Decision:* QC runs on the filtered data with its NaNs. The imputed data is used only for correction methods that need complete data. On the kNN side, hiding 2,000 known values and imputing them gave a correlation of 0.986 with the true values.

## 6. Design and structure

**6.1 Different data versions were hard to keep consistent.**
QC had to run on several versions (transformed, filtered, imputed, corrected), and it wasn't clear how to keep data and metadata aligned or compare the right versions.
*Fix:* `ProteinData` now stores one matrix per stage, and `get_stage` always returns the matching metadata. The "before correction" comparison uses the stage just before correction, so before and after cover exactly the same proteins and samples.

**6.2 Pipeline vs experimenting.**
I wanted both a full pipeline and the ability to test things one function at a time.
*Decision:* All logic stays in functions that take plain DataFrames. Pipeline wrappers will only read a stage, call the same function and store the result, so the two ways of using the package can't give different answers.

**6.3 Removed things that didn't work:**
- The silhouette score used a PCA that changed nothing (PCA without reducing dimensions is a rotation, and silhouette distance doesn't change under rotation) and filled NaN with 0, which on raw data made missingness dominate. Removed for now.
- The "batch effect detected: True/False" flag based on a single cutoff was removed, since one threshold can't decide this.

## 7. Development environment

- `pip` failed with "Unable to create process" inside the conda environment, most likely because of the space in my Windows user path. Using `python -m pip` and `python -m pytest` avoids it.
- In Jupyter, edited modules weren't reloaded, which produced errors like "unexpected keyword argument 'n_testable'" even though the code was correct. Restarting the kernel or using `%autoreload 2` fixes it.
- A test failed with "assignment destination is read-only" because newer pandas returns a read-only array from `.to_numpy()`. Fixed with `to_numpy(copy=True)`.

## 8. Next steps

1. Rewrite the pipeline wrappers using `get_stage` / `set_stage`, and remove the outdated `qc_metrics.py`.
2. Fix the `config.py` trailing-comma bug and remove unused settings.
3. Batch correction: test whether InMoose's ComBat accepts missing values. Planned options are per-plex median centering (handles missing values directly) and ComBat on imputed data with the imputed cells masked out before QC.
4. Compare QC before and after correction.
