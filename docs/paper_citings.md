
## metwarebio.com/proteomics-quality-control-reproducible-data/

## pmc.ncbi.nlm.nih.gov/articles/PMC13106926/  Practical Impact of Imputation and Batch‐Effect Correction for Proteomics/Peptidomics Differential‐Abundance Analysis (Charis Gonidaki, Agnieszka Latosinska, Antonia Vlahou, Rafael Stroggilos, Harald Mischak)

"2.4.1. Imputation Methods
Gaussian imputation replaces missing values by random numbers drawn from a normal distribution that is shifted toward lower intensities and is narrower than the empirical distribution. Following the widely accepted Perseus‐style characteristics, we sampled imputed values from a distribution with mean μ − 1.8σ and standard deviation 0.3σ, where μ and σ are the mean and standard deviation of the observed log‐intensities within each sample [20]. It is a method which assumes that missing peptide intensities correspond to low‐abundance signals just below the detection limit.

Half limit of detection (1/2 LOD) replacement sets all missing values to half the minimum observed intensity in the dataset. This approach assumes that missingness arises because true peptide abundances fall below the instrument's detection threshold [21].

K‐Nearest Neighbors (KNN) replaces a missing value by averaging the values of the K most similar features, where “similarity” is defined via a chosen distance metric [22]. In this study, we applied kNN() function from VIM package [23] to our peptide‐intensity matrix, using K = 5 and package's default Gower‐type distance (which scales numeric variables to [0,1] before computing Euclidean distance).

2.4.2. Batch Correction Methods
Three batch‐correction variations were applied to the log2‐transformed imputed datasets:

ComBat applies either parametric or non‐parametric empirical Bayes frameworks to adjust data for batch effects [24]. For this analysis, we applied the parametric empirical Bayes framework (package default) using the ComBat() function from sva package [25] on our log2‐transformed imputed peptide matrix. No biological covariates were included, so all between‐batch variations were subject to removal.

ComBat using the CKD covariate extends the above method by supplying a design matrix with the disease‐status covariate [11]. In this way, ComBat preserves differences associated with CKD status while still removing unwanted batch effects, ensuring that true case–control signal is not inadvertently washed out.

Mutual Nearest Neighbors (MNN) identifies mutual nearest neighbors between datasets to align their shared structure [26]. These neighbors are pairs of data points (samples) from different batches that are closest to each other in the high‐dimensional space, representing similar biological states. For this, we applied the function batchCorrect() from batchelor package [27] to our log2‐transformed imputed peptide matrix with PARAM = ClassicMnnParam(k = 15), and cos.norm.out = FALSE.

"

## Gaussian processes for missing value imputation (Bahram Jafrasteh, Daniel Hernández-Lobato, Simón Pedro Lubián-López, Isabel Benavente-Fernández)


