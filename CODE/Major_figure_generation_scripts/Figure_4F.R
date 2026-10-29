suppressPackageStartupMessages({
  library(data.table)
  library(scRNAtoolVis)
  library(ggplot2)
})

source_csv <- "Figure_4F_source_data.csv"
output_pdf <- "Figure_4F_volcano_reproduced.pdf"

dt <- fread(source_csv)
required_columns <- c(
  "cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "highlight"
)
stopifnot(all(required_columns %in% names(dt)))

iron_genes_B <- dt[highlight == TRUE, unique(gene)]

pal_sample <- c(
  PDX98A_arc = "#1F78B4",
  PDX98B_arc = "#E31A1C",
  PDX98C_arc_rerun = "#33A02C"
)

set.seed(20260219)

p <- jjVolcano(
  diffData = dt,
  log2FC.cutoff = 0.25,
  size = 3.5,
  fontface = "italic",
  tile.col = pal_sample,
  topGeneN = 0,
  myMarkers = iron_genes_B
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 8,
  height = 6.7,
  units = "in",
  device = "pdf"
)
