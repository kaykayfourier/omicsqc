# omicsQC

A quality control and batch-correction pipeline for bulk proteomics data, built around a real, messy dataset: ROSMAP/Banner CPTAC (Synapse: syn52854669).

Proteomics data from FragPipe tends to come with batch effects baked in different runs, different days, same underlying biology drowned out by technical noise. This package was built to find that noise, prove it's really there and correct for it without accidentally erasing the biology you care about.

## What it actually does

- Takes a FragPipe protein abundance CSV plus a metadata/traits CSV and gets them into shape
- Flags and filters out low-quality samples and proteins
- Normalizes with log2 + median scaling
- Imputes missing values in a batch-aware way — MinProb per protein (not per sample), with per-batch completeness checks for proteins missing entirely within a batch, plus a KNN-within-batch option
- Corrects for batch effects with ComBat, while explicitly protecting biological signal (diagnosis is passed in as a covariate so it doesn't get regressed out along with the batch noise)
- Validates that correction actually worked, with before/after ANOVA F-stat checks. This is how we confirmed the batch effects were real and pervasive in the first place.
- Exports clean CSVs and a markdown audit log of what happened in each run
- Comes with a small Streamlit dashboard to actually look at the results instead of squinting at a terminal

## What it doesn't do yet

Being upfront about this: there's no MissForest imputation, no GPU acceleration, no PVCA or SVA, and no multi-omics support. These are reasonable next steps, they're just not built **yet**.

## Who this is for

It was built for FragPipe-style bulk proteomics data with a known batch structure, like ROSMAP/Banner CPTAC. If your data looks like that, this should save you real time. If it doesn't, a more general-purpose tool is probably the better fit — and that's fine.

## Status

What's listed under "what it actually does" is working and validated. Batch correction layer is yet to be developed since more general-based tools like omicsGMF exist.
