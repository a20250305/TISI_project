# TISI_project

## Overview

This repository contains the core scripts and example datasets used for the development and evaluation of the Therapy-Induced Senescence Index (TISI), a transcriptome-wide framework for quantifying therapy-induced senescence across cancer types.

---

## Repository contents

### Core scripts

- `Code_Step1_Cancer_type_specific_genes.R`  
  Identification and preprocessing of cancer type-specific genes used for TISI model construction.

- `Code_Step2_OOD_and_TISI.R`  
  TISI model training, out-of-distribution (OOD) evaluation, scoring, and performance assessment.

---

### Example datasets

- `Raw_counts_test.csv`  
  Example raw count matrix used for reproducing TISI model training, testing, and OOD evaluation workflows.

- `14 datasets metadata.xlsx`  
  Metadata annotation file corresponding to the datasets used for TISI model development and evaluation.

- `FPKM_merged_MED12KO.txt`  
  Example processed expression matrix that can be directly used for TISI scoring and downstream biological validation analyses.

---

## System requirements

### Operating systems tested

- macOS
- Rocky Linux 9

### Software requirements

- R (v4.5.1 / v4.5.3)

### Major R packages

- Seurat (v5.3.0)
- harmony
- data.table
- dplyr
- ggplot2
- patchwork

Additional package dependencies are included within the scripts.

---

## Instructions for use

Users may apply the TISI framework to their own transcriptomic datasets by replacing the example expression matrices with user-generated bulk or single-cell RNA-seq data.

Input matrices should contain genes as rows and samples/cells as columns.

The `FPKM_merged_MED12KO.txt` file can be directly used as an example input for TISI scoring.

Detailed downstream analyses and figure-generation scripts are available upon reasonable request.

---

## Data availability

Raw sequencing datasets generated in this study are being deposited in GEO.

