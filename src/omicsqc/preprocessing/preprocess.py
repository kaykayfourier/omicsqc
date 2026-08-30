from omicsqc.preprocessing import impute, transform
from omicsqc.models import ProteinData
from omicsqc.config import QCconfigs

def preprocess_pipeline(protdata: ProteinData, protein_missingness_threshold_pct,sample_missingness_threshold_pct,imputation: str,batch_col: str | None = None,):

    protdata = transform.transform_data(protdata, protein_missingness_threshold_pct,sample_missingness_threshold_pct, batch_col)

    if imputation == "MinProb":
        protdata.imputed_data = impute.MinProbImpute(protdata.transformed_data, protdata.filtered_metadata, batch_col, QCconfigs.MinProb_lmbda, QCconfigs.MinProb_delta)
    else:
        print("Select appropriate imputation method. eg. [MinProb, kNN_w_batch, mean, median]")
    protdata.logs.append(f"Applied {imputation} imputation technique on the protein expression data.")
    return protdata
