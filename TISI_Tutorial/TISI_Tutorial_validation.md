# TISI tutorial validation report

Date: 2026-10-24

## Added

* TISI_Tutorial.Rmd
* TISI_Tutorial.html
* README.md
* example_data/GSE189764_CTL_vs_TIS_expression.csv.gz
* example_data/GSE189764_CTL_vs_TIS_metadata.csv
* example_data/GSE189764_CTL_vs_TIS_expected_TISI.csv
* model/Step6_SenOCLR_L2_10_final_dot.rdata
* results/GSE189764_TISI_scores.csv
* results/sessionInfo.txt
* output_figures/GSE189764_CTL_vs_TIS_TISI.pdf

The model file is byte-for-byte identical to the author-provided
Step6_SenOCLR_L2_10_final_dot.rdata. The expression example contains 38,244
genes and six samples. No simulated values were created.

## Sample verification

The example uses C1-C3 as CTL and gT1-gT3 as the TIS persister condition. The
metadata records GSM5706915-GSM5706917
and GSM5706921-GSM5706923, respectively, together with the original local
expression-file provenance.

## Verification performed

* All 24 executable R Markdown chunks ran successfully with knitr.
* The calculated scores used 22,345 mapped expression genes and 15,038 of the
  15,729 model genes.
* The maximum absolute difference between recalculated and checked-in
  TISI_scaled values was 4.440892e-16; the regression tolerance is 1e-10.
* The sample-level result CSV and one-page PDF figure were generated.
* An HTML tutorial was generated with knitr::knit2html().
* Standard rmarkdown::render() was attempted but could not run because pandoc
  is not installed in the current environment. The generated HTML is therefore
  a fallback render; installing pandoc allows the standard render command in
  README.md to be used.

## Scientific note

Figure 2F is the GSE189764 RNA-seq TISI analysis and was generated after
normalization and min-max scaling of all 12 GSE189764 samples, whereas this
public worked example intentionally ships only the required 3+3 samples.
Figure 2E is a separate C12-FDG flow-cytometry/MFI panel. The tutorial
documents this dataset-dependent scaling difference and does not claim that
the six-sample scaled values are identical to the Figure 2F values.

This tutorial addresses the worked-usage/vignette request. It does not by
itself provide the separate downstream-analysis and complete manuscript
figure-generation collection requested by the reviewer.
