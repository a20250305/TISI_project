# TISI: Therapy-Induced Senescence Index

**Code, example data, model weights, and selected analysis/figure-generation scripts accompanying the revised manuscript _“A Transcriptome-wide Therapy-Induced Senescence Index Reveals Iron-Associated Cellular Niches Underlying Chemoresistance.”_**

## Overview

The Therapy-Induced Senescence Index (TISI) is a transcriptome-wide, OCLR-derived score developed to characterize transcriptional states associated with therapy-induced senescence (TIS) across cancer contexts. This repository provides a **validated worked example for applying a pretrained TISI model**, model-development scripts, selected benchmark and figure-generation scripts, and supporting data files.

**Start here:** [TISI worked tutorial](TISI_Tutorial/README.md). The tutorial is for **inference with a pretrained model**, not retraining or complete reproduction of every analysis in the manuscript.

## Repository map

| Location | Purpose | Reproducibility status |
| --- | --- | --- |
| [`TISI_Tutorial/`](TISI_Tutorial/) | A six-sample GSE189764 bulk RNA-seq worked example; pretrained model, expression data, expected scores, and output files | Core calculation and score regression check documented as validated; standard R Markdown rendering requires Pandoc |
| [`TISI_Tutorial/model/Step6_SenOCLR_L2_10_final_dot.rdata`](TISI_Tutorial/model/Step6_SenOCLR_L2_10_final_dot.rdata) | Pretrained TISI model weights (OCLR, L2 = 10, dot-product scoring) | Provided for inference |
| [`CODE/Code_Step1_Cancer_type_specific_genes.R`](CODE/Code_Step1_Cancer_type_specific_genes.R) | Identify cancer-type-associated genes and explore lineage-gene removal | Research workflow; requires preprocessing inputs/intermediates not configured for immediate execution from a fresh clone |
| [`CODE/Code_Step2_OOD_and_TISI.R`](CODE/Code_Step2_OOD_and_TISI.R) | Cancer-type-held-out OOD assessment, OCLR model selection, and scoring | Research workflow; requires an intermediate filtered expression object and path configuration |
| [`CODE/Major_figure_generation_scripts/`](CODE/Major_figure_generation_scripts/) | Selected analytical and plotting scripts for manuscript figures | Mixed: some plotting scripts can use exported source data after path/filename configuration; some require additional data, model objects, or specialized packages |
| [`CODE/TISI_human_Xenium_plot_only.Rmd`](CODE/TISI_human_Xenium_plot_only.Rmd) and [`CODE/TISI_mouse_Xenium_plots_only.ipynb`](CODE/TISI_mouse_Xenium_plots_only.ipynb) | Targeted spatial-transcriptomics visualization | Plotting-oriented materials, not a stand-alone raw-data processing pipeline |
| [`Data/Figure_generation_source_data/`](Data/Figure_generation_source_data/) | Source-data files for selected figure panels, including CSV and compressed CSV | Use as supplied panel-level data; not a substitute for the complete upstream analysis inputs |
| [`Data/14 datasets metadata.xlsx`](Data/14%20datasets%20metadata.xlsx) | Development cohort metadata | Supporting metadata |
| [`Data/Raw_counts_test.csv`](Data/Raw_counts_test.csv) | Example count matrix | Research example; not a substitute for every intermediate required by the training scripts |
| [`Data/Preprocessed raw counts matrix of 14 datasets_part1.csv`](Data/Preprocessed%20raw%20counts%20matrix%20of%2014%20datasets_part1.csv), [part2](Data/Preprocessed%20raw%20counts%20matrix%20of%2014%20datasets_part2.csv), [part3](Data/Preprocessed%20raw%20counts%20matrix%20of%2014%20datasets_part3.csv) | Split preprocessed count-matrix files for the 14 development datasets | Combine and inspect according to the matrix headers/metadata before use; the current scripts do not automatically reconstruct all necessary inputs |
| [`Data/FPKM_merged_MED12KO.txt`](Data/FPKM_merged_MED12KO.txt) | Processed expression matrix for the MED12-KO use case | Example downstream input; preprocessing must be consistent with the intended scoring workflow |

## Quick start: applying the pretrained TISI model

1. Clone/download this repository, and open a terminal in `TISI_Tutorial/`.
2. Install R and the tutorial's dependencies. **Pandoc** is required for `rmarkdown::render()` (it is generally available with RStudio).
3. Run the R Markdown tutorial.

```bash
git clone https://github.com/a20250305/TISI_project.git
cd TISI_project/TISI_Tutorial
```

```r
install.packages(c("data.table", "ggplot2", "rmarkdown", "knitr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
rmarkdown::render("TISI_Tutorial.Rmd")
```

The tutorial loads six GSE189764 samples (three controls and three TIS-persister samples), matches input genes to the pretrained model, calculates scores, compares the computed scaled values with the checked-in expected-results file, and writes:

- `TISI_Tutorial/results/GSE189764_TISI_scores.csv`
- `TISI_Tutorial/output_figures/GSE189764_CTL_vs_TIS_TISI.pdf`

See the [tutorial README](TISI_Tutorial/README.md), [executable R Markdown](TISI_Tutorial/TISI_Tutorial.Rmd), and [tutorial validation notes](TISI_Tutorial/TISI_Tutorial_validation.md) for the exact dataset, tested environment, and validation results. The validation notes report calculation-level checks; they also note that a standard R Markdown render was not tested successfully in the original environment because Pandoc was absent.

## Interpreting TISI scores

- The pretrained model uses **OCLR weights with L2 = 10** and **dot-product scoring**.
- Input data must be a **gene × sample** matrix (or gene × cell matrix for properly validated per-cell scoring), with gene identifiers that can be mapped to the model's Ensembl IDs.
- The worked tutorial uses a specific normalized expression representation and applies log transformation, column-wise scaling, gene-wise centering, and the weighted dot product. **Raw counts, TPM/FPKM, and already log-transformed values are not interchangeable.** Follow the documented preprocessing for the intended workflow rather than substituting arbitrary expression scales.
- `TISI_raw` is the model-based score; `TISI_scaled` is a **within-analysis min–max transformation to [0, 1]**. Because preprocessing and scaling can depend on the analyzed sample set, values computed or min–max scaled independently across datasets are **not automatically on a shared absolute scale**.
- A higher TISI score indicates stronger similarity to the model's TIS-associated transcriptional state. It is **not a probability, a diagnostic classification, or a validated universal clinical cutoff**.
- Always report the number/fraction of model genes matched and the preprocessing procedure when transferring the model to another platform.

The tutorial is demonstrated on bulk RNA-seq. Its optional single-cell template is illustrative and **should not be treated as a validated replacement for the manuscript's full single-cell or spatial preprocessing workflow** without checking normalization, input layer, gene coverage, and sample-level aggregation.

## Model development and out-of-distribution evaluation

The manuscript describes a pan-cancer development set of **258 samples from 40 cell lines, 12 cancer types, and 14 source datasets** (152 senescent and 106 control samples). The model-development code is contained in:

1. [`CODE/Code_Step1_Cancer_type_specific_genes.R`](CODE/Code_Step1_Cancer_type_specific_genes.R) — marker discovery, ARI sweep, lineage-related gene filtering.
2. [`CODE/Code_Step2_OOD_and_TISI.R`](CODE/Code_Step2_OOD_and_TISI.R) — OOD cross-validation, penalty and scoring-method exploration, and fitted-model/scoring outputs.

The Step 2 script constructs **three shuffled repetitions of ten cancer-type-held-out folds (30 test partitions)**. These are not 30 statistically independent source cohorts. It evaluates six L2 settings and multiple scoring approaches. The configuration reported in the revised manuscript uses **L2 = 10 with dot-product scoring**.

**Execution caveat:** The two scripts are released as research workflows and are **not currently an end-to-end turnkey training pipeline from a fresh clone**. For example, Step 1 refers to `Counts_Global.csv`, `raw_full.rdata`, and `DEGs_CancerType_list.rds`; Step 2 expects an intermediate filtered expression object at `./data/exprData_dropCancerType_Top50.rdata`. These inputs/locations must be reconciled with the provided data files and the scripts' output paths. Verify the preprocessing and model-selection steps before treating a newly trained model as numerically identical to the supplied pretrained model.

## Benchmarking and figure-generation scripts

- Fig. 2A–B benchmarking: [`Figure_2A_benchmark.R`](CODE/Major_figure_generation_scripts/Figure_2A_benchmark.R) and [`Figure_2B_benchmark.R`](CODE/Major_figure_generation_scripts/Figure_2B_benchmark.R). The manuscript compares TISI against **31 other methods**, with **8 bulk transcriptomic datasets** (RNA-seq and microarrays) and **9 single-cell comparisons from 6 GEO series**.
- Additional scripts for selected panels of Figs. 1–6 and Fig. S1 are in [`CODE/Major_figure_generation_scripts/`](CODE/Major_figure_generation_scripts/).
- Related panel-level files are in [`Data/Figure_generation_source_data/`](Data/Figure_generation_source_data/); consult actual filenames and column headers before running a plotting script.
- Complex spatial analyses, including SpatialCellChat, require upstream expression, location/deconvolution objects, and (where applicable) a specific ligand–receptor database. The source scripts expose these inputs but **do not package all of them in this repository**.

**Execution caveat:** Several figure scripts contain project-specific relative paths, placeholder paths, or filenames that differ from the exported CSV names. For example, `Figure_4J.R` expects `Figure_4J_panel1_source_data.csv`, while the public source-data folder contains `Figure_4J_source_data.csv`. Adapt file paths and confirm that the provided source-data columns correspond to the plotting code before execution. Availability of a panel-level CSV should not be interpreted as complete independent reconstruction from raw sequencing data.

## Software and environment

- **Pretrained-model tutorial, documented validation:** R 4.5.1 on macOS; package versions are recorded in [`TISI_Tutorial/results/sessionInfo.txt`](TISI_Tutorial/results/sessionInfo.txt).
- **Original manuscript training environment:** R 4.0.5 is reported in the Methods. Some current public scripts use modern package APIs (including Seurat v5), so the original manuscript environment and the public tutorial environment are **not interchangeable**.
- **Core development dependencies:** `dplyr`, `DESeq2`, `Seurat`, `mclust`, `ggplot2`, `gelnet`, `GSVA`, `pbapply`, `pROC`, and packages required by the relevant scripts.
- **Figure analysis/plotting:** additional libraries such as `data.table`, `ComplexHeatmap`, `circlize`, `patchwork`, `SpatialCellChat`, and others may be needed, depending on the panel.

Package versions are not presently pinned across every workflow. Record `sessionInfo()` and any script-specific inputs for each analysis. The tutorial supplies a concrete test case; a fully automated clean-environment reproduction test for every figure is not asserted here.

## Data provenance and availability

Public source accessions and cohort definitions should be taken from the manuscript's supplementary dataset tables. The files under `Data/` provide selected matrices, metadata, and plotted source data; the repository **does not contain every original raw sequencing run, upstream intermediate object, or specialized spatial database** needed to regenerate all manuscript figures. Consult the manuscript's Data Availability statement for the status/accessions of newly generated sequencing data, and add specific accession identifiers or archived release links here once confirmed.

## Reproducibility and support

The most direct entry point for a new user is the [validated six-sample inference tutorial](TISI_Tutorial/README.md). For reproducing a manuscript panel, begin with the matching script in [`CODE/Major_figure_generation_scripts/`](CODE/Major_figure_generation_scripts/) and its exported source data, then inspect all required external inputs. For questions or discrepancies, please open a GitHub issue with the script, input schema, R/package versions, and full error message.

*This README documents the repository as currently distributed; it deliberately distinguishes a tested inference demonstration from research workflows requiring further input preparation and environment configuration.*
