# Proposed Metrics to be used to assess Batch Effect
1. Per-protein ANOVA (or Kruskal-Wallis)
    Answers the question "Is there evidence that batch affects the data?"
2. Mean R² of batch + Median R² + Distribution of R² + PVCA (if implemented)
    Answers the question "How much variation is explained by batch?"
3. Compare R²(batch) vs R²(disease) or use variance partitioning.

4. Silhouette score (batch labels)

5. Chi-square test/Logistic regression to determine missingness to be batch dependent 

Run these fundamental metrics before and after batch correction


## ANOVA Per Protein statistical test
- ANOVA test assumes a normal distribution in your dataset. Additionally, the ANOVA statistical tests is fairly robust against normality deviations given that the sample sizes are sufficiently large. 
- ANOVA assumes that the variances of the populations that the samples come from are equal.
- It is recommended that the user apply a log2 Transform and Normalization transformation on the dataset before calling compute_batch_metrics() as its default method is ANOVA.

## Kruskal-Wallis Per Protein test
- Kruska-Wallis is a non-parametric test that serves as an appropriate alternative to the ANOVA test given that the assumptions of ANOVA test are violated (absence of normal distribution & redundant values skewing the distribution).
- It employs a ranking method and uses a tie correction called Mid Ranks to effectively compute ranks.
- User must explicitly call method = "kruskal_wallis" to perform this test given that their dataset violates ANOVA assumptions to get reliable results.


## Variance metrics
Functions:
1 -> variance_explained() {
    return mean_r2 + med_r2 + std_r2 + q25 + q75 + max_r2 + min_r2
} 
2 -> compare_variance_explained(X1_variance_explained, X2_variance_explained) {}
later additions:
3 -> PVCA()

4 -> variance_partitioning()

