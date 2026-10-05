## Dataset assumptions:
- Module takes in only 2 datasets: an Expression dataset and Metadata dataset.
- Structure of Expression dataset is protein_id x sample_id
- Structure of Metadata dataset is sample_id x metafeatures
- The defualt naming convention follows the TMT ROSMAP and BNNR, where sample_id is: <Cohort>_<Batch>.<Channel> for example `rmap_b32.129N` and `bnnr_b14.130N`
- The specific Batch column required by the pipeline must be passed by the user with the knowledge that it adheres to the same naming convention of <Cohort_Batch/Plex> like `rmap_b32` to indicate the right technical batch column. Failure to do so will throw format convention error and will result in an interrupted data ingestion. 
- In case the Expression dataset has no protein id index, the ingestor will fallback to the first column as the protein_id column.

## Batch Metrics
### Batch vs. cohort confounding

In ROSMAP/Banner TMT data every batch (plex) belongs to exactly one cohort (rmap_* or bnnr_*). The batch test compares plexes as given, so differences between cohorts (biology, sample handling, or processing) are counted as "batch effect". A large F or eta² therefore means "plexes differ", not necessarily "the instrument run caused it".

To separate the two, run the test within each cohort (subset the samples to one cohort and call the test again), and compare the eta² distributions. Signal that persists inside a cohort is batch; signal that appears only when cohorts are mixed is partly cohort.

min_batches (default 3) excludes proteins seen in too few plexes; excluded proteins are marked low_coverage in the result and listed in the log.

**NOTE TO THE USER:** min_batches defaults to 3. Proteins seen in only a few batches can still produce very large F values from few samples; raise it (e.g. to 25% of the number of batches) for stricter filtering.

## Missingness Analysis
- Batch as a covariate is true by construction. With partial_fraction near 0, every sample in a plex shares one missingness value, so the batch row (effect size 0.9998 on your data) is expected and carries no information.
- Sample-level p-values are inflated. The 8-9 samples in a plex are one outcome, so the real sample size for missingness is the number of plexes (67), not 558. For other covariates, read effect_size, not the p-value.
- Missingness is plex-level on this dataset (0.0% partly-missing cells), so the intensity result is about detection per plex.
- MAR vs MNAR is not identifiable from the observed data. The verdict only reports intensity dependence.

## Variance explained (variance layer)

For each protein, this computes the adjusted R² of a categorical factor (batch, diagnosis, sex, ...), so batch can be compared with the biology it may hide. Run it on transformed_data and again on corrected_data. Input must be log2-transformed and median-centered.

Why adjusted R². Plain R² rises with the number of groups even for pure noise (about (k-1)/(n-1)), so a 67-level batch would look larger than a 2-level diagnosis. Adjusted R² = 1 - (SS_within/(n-k)) / (SS_total/(n-1)) corrects for this. Zero means the factor explains no more than chance, and values can be slightly negative.

The cogdx issue. cogdx is recorded only for the rmap cohort (all bnnr samples are missing), so it was being measured on different samples than the other factors. Proteins seen in only a few plexes also inherit large plex offsets (about 5 log2 units), which produced cogdx R² of 0.77-0.93 on a few proteins. Two safeguards fix this:

complete_cases=True (default): every factor is computed on the same samples, those with a value for all factors. The sample count is logged and reported as n_samples.

min_obs=30 (default): proteins observed in fewer samples are not tested. This is a heuristic and can be tuned.

Limits:

 Values are per factor and factors overlap (batch is nested in cohort), so they should not be added up. Before correction, the R² of a biological factor can include batch signal if the factor is unevenly spread across plexes. With complete cases including cogdx, batch is measured on rmap only (360 samples, 45 plexes). After correction, expect batch R² to fall and biological R² to rise (the same effect becomes a larger share of the remaining variance).

## Filtering: choosing the default threshold

The old filter dropped any protein that was fully missing in at least one batch. Missingness in ROSMAP/Banner TMT data is plex-level (0% of protein-plex cells are partly missing), so this removed nearly every protein with any missing value and left almost nothing to impute. That rule was removed. Proteins are now filtered only on overall missingness, after samples.

We tested the default (proteins ≤ 80% missing, samples ≤ 50% missing) on the full dataset (558 samples, 11,958 proteins):

- All 558 samples and all 67 batches were kept (8-9 samples each).
- 10,448 proteins (87%) were kept. Missingness fell from 22.3% to 12.3%, so the imputation step still has work to do.
- The 1,510 dropped proteins had a median of 92.7% missing, i.e. seen in about 5 of 67 plexes.
- Low-abundance proteins are thinned but not removed: 40% of the lowest intensity decile and 70% of the next one were kept, and over 99% of proteins above that.

A threshold sweep showed no clear elbow:

| Max protein missing | Proteins kept | Missing left |
|---|---|---|
| 20% | 8,152 | 1.6% |
| 50% | 9,304 | 5.7% |
| 80% | 10,448 | 12.3% |
| 95% | 11,374 | 18.5% |

So 80% was chosen as a lenient default that keeps most low-abundance proteins for imputation while dropping the sparsest ones. On this dataset it means "observed in at least ~14 of 67 plexes". It is a parameter: 50% is the stricter option if less imputation is preferred.

## Imputation
The omicsqc.preprocessing.impute imports 3 functions: impute_knn(kNN: for values missing at random.), impute_minprob (MinProb (Lazar et al., imputeLCMD): for left-censored values) and impute_mixed (groups random samples and left-censored samples; applies kNN and MinProb respectively)

**NOTE TO USER:** It has been observed that the eta^2 of Imputed data slightly reduces due to MinProb imputation, credited to MinProb diluting the batch effect for censored proteins. It is then recommended not to run QC metrics on filtered, log2 transformed dataset and not imputed dataset for more reliable QC metrics. The imputation serves as a prerequisite to batch correction and not pre-correction QC analysis.