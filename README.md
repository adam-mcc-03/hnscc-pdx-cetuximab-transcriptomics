# Transcriptomic Analysis of Cetuximab Response in HNSCC PDX Models

This repository contains the principal R scripts supporting the MSc thesis:

“Transcriptomic Analysis of Cetuximab Response in Head and Neck Squamous Cell Carcinoma Patient-Derived Xenograft Models.”

Adam McClatchey

MSc Bioinformatics and Computational Genomics

Queen’s University Belfast, 2026

## Scripts
01_Cohort_Response_QC.R
Cohort characterisation, response classifications, trajectory analysis and baseline quality control.

02_Baseline_DE_Pathway.R
Baseline differential expression and pathway enrichment analyses, including cross-classification DEG and pathway recurrence.

03_TF_Activity_Master_Regulators.R
DoRothEA/VIPER transcription-factor activity analysis, TF recurrence and TF–DEG mapping for candidate master regulators.

04_Longitudinal_Analysis.R
Mixed-effects differential expression and pathway analyses across untreated, short-treatment and endpoint samples.

## Analytical structure

The complete baseline differential-expression, pathway-enrichment and transcription-factor pipelines are presented for the representative average 3–6-week responder versus non-responder comparison. The same analytical workflows were repeated for the remaining response classifications. Their saved outputs were subsequently integrated in the cross-classification analyses.

The scripts use relative file paths and should be run from the main project directory containing the required input and previously generated result folders. Required R packages are listed at the beginning of each script.

## Data availability

The input RNA-sequencing data, metadata and model-level result files are not included in this repository. These data were provided for the MSc research project and are retained by the research group; public redistribution has not been authorised. The repository is solely for the purpose of proof of code. 

This repository is provided to document the analytical methods and code used in the thesis rather than as a self-contained data and software package.
