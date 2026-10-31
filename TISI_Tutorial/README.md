# TISI tutorial: GSE189764 worked example

This folder is a self-contained worked example for applying the pretrained
Therapy-Induced Senescence Index (TISI) model to bulk RNA-seq data.

## Quick start

From this directory:

~~~r
install.packages(c("data.table", "ggplot2", "rmarkdown", "knitr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
rmarkdown::render("TISI_Tutorial.Rmd")
~~~

The tutorial loads the six-sample GSE189764 test dataset, calculates the
pretrained L2=10 dot-product TISI score, checks the result against the expected
CSV, exports sample-level scores, and creates a PDF plot.

## Files

* [TISI_Tutorial.Rmd](TISI_Tutorial.Rmd): executable step-by-step walk-through.
* [example_data/GSE189764_CTL_vs_TIS_expression.csv.gz](example_data/GSE189764_CTL_vs_TIS_expression.csv.gz): all-gene expression matrix, genes x 6 samples.
* [example_data/GSE189764_CTL_vs_TIS_metadata.csv](example_data/GSE189764_CTL_vs_TIS_metadata.csv): sample-to-condition mapping and GEO accessions.
* [example_data/GSE189764_CTL_vs_TIS_expected_TISI.csv](example_data/GSE189764_CTL_vs_TIS_expected_TISI.csv): expected raw and min-max-scaled scores.
* [model/Step6_SenOCLR_L2_10_final_dot.rdata](model/Step6_SenOCLR_L2_10_final_dot.rdata): pretrained named Ensembl gene weights used for inference.
* `output_figures/`: formal PDF figure outputs from the tutorial.
* `knitr_assets/`: PNG assets used by the rendered HTML/Markdown tutorial.

The example uses C1-C3 as CTL and gT1-gT3 as the TIS persister condition. The
source series is GSE189764; the six GEO
sample accessions are recorded in the metadata CSV.

## Input requirements

The required orientation is genes x samples (or genes x cells for the optional
scRNA-seq template). The original RNA-seq pipeline accepts the normalized
expression representation used for this example and applies log2(x + 1e-5),
column-wise scaling, row-centering, model-weight matching, and a dot product.
The tutorial documents duplicate-gene averaging, Ensembl-version cleanup,
missing-gene handling, model coverage, raw scores, and min-max scaling.

## Interpretation

TISI_raw is the weighted dot-product score. TISI_scaled is a dataset-level
0-to-1 min-max transformation. Higher values indicate a more TISI-like
transcriptional state relative to the samples used for that calculation. The
score is not a probability and is not a universal clinical cutoff.

The original Figure 2F is the GSE189764 RNA-seq TISI analysis and used all 12
GSE189764 samples during normalization and min-max scaling before displaying
the six samples. Figure 2E is a separate C12-FDG flow-cytometry/MFI panel. The
tutorial intentionally ships only the required 3+3 samples, so its six-sample
scaled values differ slightly from the Figure 2F values. This reflects
dataset-dependent scaling, not a change in the TISI model. Use the complete
12-sample matrix to reproduce Figure 2F numerically.

## Model training versus inference

This tutorial only loads the pretrained model and performs inference. It does
not retrain the OCLR model or estimate new gene weights. Downstream analyses
and manuscript figure-generation scripts are separate from this minimal usage
tutorial and should be linked separately in the public repository.
