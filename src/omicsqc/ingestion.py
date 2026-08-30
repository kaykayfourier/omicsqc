import numpy as np
import pandas as pd
from .models import ProteinData


def standardize_protein_index(data_df: pd.DataFrame) -> pd.DataFrame:
    if isinstance(data_df.index, pd.RangeIndex) or data_df.index.name is None:
        first_col = data_df.columns[0]
        data_df = data_df.set_index(first_col)
        
    data_df.index.name = "uniprot_id"
    return data_df

def ingest_data(protein_data: pd.DataFrame, meta_data: pd.DataFrame, float32_compress = True ) -> ProteinData:
    
    #Setting Index of Proteins
    data_df = standardize_protein_index(protein_data)

    #Identifying Index of samples in MetaData
    possible_id_cols = ['Unnamed: 0', 'sample_id', 'specimenID', 'sample', 'Run', 'id']
    meta_index = next((col for col in meta_data.columns if col in possible_id_cols or 
                        any(kw in col.lower() for kw in ['id', 'sample', 'specimen'])), meta_data.columns[0])


    #Memory reduction
    if float32_compress == True:
        data_df = data_df.astype(np.float32)
    #Inner join of common samples
    data_samples = set(data_df.columns)
    meta_samples = set(meta_data[meta_index])
    common_samples = list(data_samples.intersection(meta_samples))
    print(f"[Success] Found {len(common_samples)} overlapping technical samples.")
    data_filtered = data_df[common_samples].copy()
    
    #Index recheck
    meta_filtered = meta_data[meta_data[meta_index].isin(common_samples)]
    meta_filtered = meta_filtered.set_index(meta_index)
    meta_filtered = meta_filtered.reindex(data_filtered.columns)
    if meta_filtered.index.has_duplicates:
        raise ValueError(
            "Duplicate sample IDs found in metadata."
        )
    if data_filtered.columns.has_duplicates:
        raise ValueError(
            "Duplicate sample IDs found in expression matrix."
        )
    #Replacing Zeros with NaNs
    data_filtered = data_filtered.replace(0.0, np.nan)
    meta_filtered = meta_filtered.replace(0.0, np.nan)
    
    
    return ProteinData(data_filtered, meta_filtered)