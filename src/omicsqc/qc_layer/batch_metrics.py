import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from ..config import QCconfigs
from ..models import BatchMetrics

# ANOVA method Vectorized function
def anova_batch_vectorized(data_df:pd.DataFrame, meta_df : pd.DataFrame, batch_col, alpha=0.05):
    from scipy.stats import f as f_dist
    
    expr = data_df.values.astype(float)
    batch = meta_df.loc[data_df.columns, batch_col].to_numpy()
    unique_batches = np.unique(batch[~pd.isna(batch)])
    
    n_proteins, n_samples = expr.shape
    total_n = np.sum(~np.isnan(expr), axis=1)
    grand_mean = np.divide(np.nansum(expr, axis=1), total_n, out=np.full(n_proteins, np.nan), where=total_n > 0)
    
    k = len(unique_batches)
    n = np.sum(~np.isnan(expr), axis=1)
    
    ss_between = np.zeros(n_proteins)
    ss_within = np.zeros(n_proteins)
    
    for b in unique_batches:
        mask = batch == b
        group = expr[:, mask]
        group_n = np.sum(~np.isnan(group), axis=1)
        group_mean = np.divide(np.nansum(group, axis=1), group_n, out=np.full(n_proteins, np.nan), where=group_n > 0)
        valid = group_n > 0
        ss_between[valid] += (group_n[valid] * 
                             (group_mean[valid] - grand_mean.squeeze()[valid]) ** 2)
        ss_within += np.nansum((group - group_mean[:, None]) ** 2, axis=1)
    
    df_between = k - 1
    df_within = n - k
    
    with np.errstate(divide='ignore', invalid='ignore'):
        ms_between = np.where(df_between > 0, ss_between / df_between, np.nan)
        ms_within = np.where(df_within > 0, ss_within / df_within, np.nan)
        
        f_stat = np.where(ms_within > 0, ms_between / ms_within, np.nan)
    
    p_values = f_dist.sf(f_stat, df_between, df_within)
    
    anova_df = pd.DataFrame({
        "protein": data_df.index,
        "F": f_stat,
        "p_value": p_values,
        "significant": p_values < (alpha / n_proteins)
    }).set_index("protein")
    
    anova_sorted = anova_df.sort_values(by="F", ascending=False)
    anova_sorted["Flag"] = (
        (anova_sorted["F"] > anova_sorted["F"].mean()) & 
        (anova_sorted["p_value"] < 0.05)
    )
    return anova_sorted
    
# Kruskal-Wallis Test
from scipy.stats import kruskal

def kruskal_batch(data_df:pd.DataFrame, meta_df : pd.DataFrame, batch_col: str, alpha: float = QCconfigs.alpha,):
    
    print("Running Kruskal Wallis test per protein...")
    batch = meta_df.loc[data_df.columns, batch_col].to_numpy()
    unique_batches = np.unique(batch[~pd.isna(batch)])

    H = np.empty(len(data_df), dtype=float)
    P = np.empty(len(data_df), dtype=float)

    values = data_df.to_numpy(dtype=float)

    
    batch_masks = [batch == b for b in unique_batches]
    print("Looping over proteins..")
    for i in range(values.shape[0]):
        protein = values[i]

        groups = []

        for mask in batch_masks:
            g = protein[mask]
            g = g[~np.isnan(g)]

            if len(g):
                groups.append(g)

        if len(groups) < 2:
            H[i] = np.nan
            P[i] = np.nan
            continue

        try:
            H[i], P[i] = kruskal(*groups)
        except ValueError:
            H[i] = np.nan
            P[i] = np.nan

    result = pd.DataFrame(
        {
            "H": H,
            "p_value": P,
        },
        index=data_df.index,
    )

    bonf = alpha / len(result)

    result["significant"] = result["p_value"] < bonf

    #result["flag"] = result["significant"]
    print("Kruskal Wallis test successfully executed and returned.")
    return result


def sh_score(data_df:pd.DataFrame, meta_df : pd.DataFrame, batch_label: str,):
    X = data_df.copy().T.fillna(0)
    meta = meta_df.copy()
    meta = meta.loc[X.index]
    pca = PCA()
    X_pca = pca.fit_transform(X)
    labels = meta[batch_label]

    score = silhouette_score(X_pca, labels= labels)
    return score, X_pca

# Single Function to run Batch Metrics
def batch_variation_test(data_df:pd.DataFrame, meta_df : pd.DataFrame, protein_threshold_pct: np.float32, batch_effect_method: str = QCconfigs.batch_effect_method, alpha: float = QCconfigs.alpha, ):

    batch_effect_flag = False

    possible_batch_cols = ["Batch_overall","batch","Batch","Batches","MS_Batch","RunBatch"]
 
    batch_col = next(
        (c for c in meta_df.columns
         if c in possible_batch_cols or "batch" in c.lower()),
        None
    )
    #Performing Batch Effect Flagging
    if batch_effect_method == "anova":
        batch_effect_per_protein = anova_batch_vectorized(data_df, meta_df, batch_col  = batch_col, alpha = alpha)
    elif batch_effect_method == "kruskal_wallis":
        batch_effect_per_protein = kruskal_batch(data_df, meta_df, batch_col  = batch_col, alpha = alpha )
        
    else:
        print("Please choose appropriate batch effect method. [anova, kruskal_wallis]")
        return None
    
    n_proteins = batch_effect_per_protein["significant"].sum()
    n_proteins_pct = (n_proteins / len(batch_effect_per_protein)) * 100

    if n_proteins_pct >= protein_threshold_pct:
        batch_effect_flag = True

    score, data_pca_coords = sh_score(data_df, meta_df, batch_label= batch_col)
    
    batch_metrics = BatchMetrics(batch_effect_per_protein= batch_effect_per_protein,
                                 n_significant= n_proteins,
                                 n_significant_pct= n_proteins_pct,
                                batch_effect_flag= batch_effect_flag,
                                silhoutte_score= score,
                                data_pca_coords=data_pca_coords)
    
    logs = (
        f"""
    [Batch Variation Test]
    Statistical method:             {batch_effect_method} \n
    Significance level (alpha):     {alpha} \n
    Protein threshold:              {protein_threshold_pct:.2f}% \n
    Significant proteins:           {n_proteins} \n
    Significant proteins (%):       {n_proteins_pct:.2f}% \n
    Batch effect detected:          {batch_effect_flag} \n
    silhoutte score:                {score:.4f} \n
    """.strip()
    )
    return batch_metrics, logs
