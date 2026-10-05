import numpy as np
import pandas as pd
from .models import ProteinData

TMT_BATCH_RE  = r"[A-Za-z]+_b\d+"                        # rmap_b32
TMT_SAMPLE_RE = rf"(?P<batch>{TMT_BATCH_RE})\.(?P<channel>\d{{3}}[NC]?)"  # rmap_b32.129N


def standardize_protein_index(data_df: pd.DataFrame) -> pd.DataFrame:
    if isinstance(data_df.index, pd.RangeIndex) or data_df.index.name is None:
        first_col = data_df.columns[0]
        data_df = data_df.set_index(first_col)    
    
    return data_df

# TMT Naming convention check parser
def _show(s, n=3):
    return s.head(n).tolist()


def TMT_convention_validator(batch_col: pd.Series, sample_id: pd.Series):
    batch_col = pd.Series(batch_col).astype(str).reset_index(drop=True)
    sample_id = pd.Series(sample_id).astype(str).reset_index(drop=True)

    if len(batch_col) != len(sample_id):
        raise ValueError("batch column and sample ids differ in length")
    bad = batch_col[~batch_col.str.fullmatch(TMT_BATCH_RE)]
    if not bad.empty:
        raise ValueError(f"Batch column violates TMT '<cohort>_<batch>' format. eg: {_show(bad)}")

    parsed = sample_id.str.extract(f"^{TMT_SAMPLE_RE}$")
    bad = sample_id[parsed["batch"].isna()]
    if not bad.empty:
        raise ValueError(f"Sample IDs violate TMT '<cohort>_<batch>.<channel>' format, e.g. {_show(bad)}")

    mismatch = batch_col != parsed["batch"]
    if mismatch.any():
        raise ValueError(
            f"Batch column disagrees with batch parsed from sample ID, "
            f"e.g. {list(zip(_show(batch_col[mismatch]), _show(sample_id[mismatch])))}"
        )
    return True



def ingest_data(protein_data: pd.DataFrame, meta_df: pd.DataFrame, batch_col: str, sample_id_indexed = False, convention: str = "TMT", float32_compress = True ) -> ProteinData:
    
    #Setting Index of Proteins
    data_df = standardize_protein_index(protein_data)

    #Identifying Index of samples in MetaData
    if not sample_id_indexed:
        expr_cols = set(data_df.columns)
        candidates = {"<index>": meta_df.index.to_series(), **{c: meta_df[c] for c in meta_df.columns}}
        scores = {name: s.astype(str).isin(expr_cols).sum() for name, s in candidates.items()}
        meta_index_col = max(scores, key=scores.get)
        if meta_index_col != "<index>":
            meta_df = meta_df.set_index(meta_index_col)

    if batch_col not in meta_df.columns: 
        raise ValueError(f"Given batch_col : {batch_col} does not exist in metadata columns.")
    # id redundancy check
    if meta_df.index.has_duplicates:
        raise ValueError("Duplicate Metadata Sample IDs detected")
    if data_df.columns.has_duplicates:
            raise ValueError("Duplicate Expression Data Sample IDs detected")

    
    #Inner join of common samples
    meta_samples = set(meta_df.index)
    common_samples = [c for c in data_df.columns if c in meta_samples]

    if len(common_samples) == 0:
        raise ValueError("No common Sample IDs found between Expression data and Metadata")
    
    data_filtered = data_df[common_samples].copy()
    meta_filtered = meta_df.loc[common_samples]

    # less than 2 distinct batches error detection
    counts = meta_filtered[batch_col].value_counts()
    if len(counts) < 2:
        raise ValueError("Need at least 2 batches after matching samples.")
    if (counts < 2).any():
        raise ValueError(f"Batches with <2 samples after matching: {counts[counts < 2].index.tolist()}")
    # naming convention check
    if convention == "TMT":
            TMT_convention_validator(batch_col= meta_filtered[batch_col], sample_id= meta_filtered.index)

    else:
        raise ValueError(f"Unknown convention '{convention}'")

    print(f"[SUCCESS] {len(common_samples)} overlapping samples found")
    
    #Replacing Zeros with NaNs
    data_filtered = data_filtered.replace(0.0, np.nan)

    
    #Memory reduction
    if float32_compress == True:
        data_filtered = data_filtered.astype(np.float32)
    
    return ProteinData(abundance_matrix=data_filtered, metadata= meta_filtered, batch_col= batch_col)