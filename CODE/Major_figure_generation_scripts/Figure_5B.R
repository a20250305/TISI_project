suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source_csv <- "Figure_5B_source_data.csv.gz"
output_pdf <- "Figure_5B.pdf"

dt <- fread(source_csv, check.names = FALSE)
niche_cols <- paste0("Niche_", 1:13)
stopifnot(all(c("sample_label", "paper_group", "sample_spot_order", niche_cols) %in% names(dt)))
stopifnot(all(is.finite(as.matrix(dt[, ..niche_cols]))))

sample_levels <- unique(dt$sample_label)
dt[, sample_label := factor(sample_label, levels = sample_levels)]
dt[, paper_group := factor(
  paper_group,
  levels = c("Primary", "Pre-MRD", "MRD", "Recurrent")
)]
setorder(dt, sample_label, sample_spot_order)

niche_mat <- t(as.matrix(dt[, ..niche_cols]))
rownames(niche_mat) <- niche_cols
colnames(niche_mat) <- paste(dt$sample_label, dt$sample_spot_order, sep = "__")

sample_cols <- c(
  "Primary-1" = "#5B7D55", "Primary-2" = "#4B2B63", "Primary-3" = "#222222",
  "Primary-4" = "#5B6259", "Primary-5" = "#000000", "Primary-6" = "#847A7C",
  "Primary-7" = "#1D1D1D", "Primary-8" = "#91B7A2", "Primary-9" = "#E8F2D9",
  "Primary-10" = "#A6B9F1", "Pre-MRD-1" = "#4C0B50", "Pre-MRD-2" = "#AC8ACA",
  "MRD-1" = "#9C7954", "MRD-2" = "#71616D", "MRD-3" = "#3F3E70",
  "MRD-4" = "#6B6D77", "Recurrent-1" = "#E51E3D", "Recurrent-2" = "#2D4E4E"
)
group_cols <- c(
  "Primary" = "#4DBBD5", "Pre-MRD" = "#00A087",
  "MRD" = "#E64B35", "Recurrent" = "#3C5488"
)

niche_col_fun <- colorRamp2(
  c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
  c("#F6EFF7", "#E7DCEF", "#D4BFE0", "#B998C8", "#8D5DA5", "#54278F")
)

top_anno <- HeatmapAnnotation(
  Sample = dt$sample_label,
  Group = dt$paper_group,
  col = list(Sample = sample_cols, Group = group_cols),
  annotation_name_side = "left",
  annotation_height = unit.c(unit(4, "mm"), unit(4, "mm"))
)

pdf(output_pdf, width = 15, height = 5)
ht <- Heatmap(
  niche_mat,
  name = "Niche activity",
  col = niche_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9),
  top_annotation = top_anno,
  column_split = dt$sample_label,
  border = FALSE,
  use_raster = FALSE
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
