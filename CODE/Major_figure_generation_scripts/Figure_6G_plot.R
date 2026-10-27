suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source_dir <- "."
macro_csv <- file.path(source_dir, "Figure_6G_source_data_macrophages.csv")
fibro_csv <- file.path(source_dir, "Figure_6G_source_data_fibroblasts.csv")
manifest_csv <- file.path(source_dir, "Figure_6G_gene_manifest.csv")

macro_dt <- fread(macro_csv)
fibro_dt <- fread(fibro_csv)
manifest <- fread(manifest_csv)
stopifnot(all(c("spot", "sample", "TISI_scaled") %in% names(macro_dt)))
stopifnot(all(c("spot", "sample", "TISI_scaled") %in% names(fibro_dt)))
stopifnot(all(c("Module", "Gene", "GeneOrder") %in% names(manifest)))

gene_sets <- list(
  Macrophages = manifest[Module == "Macrophages"][order(GeneOrder), Gene],
  Fibroblasts = manifest[Module == "Fibroblasts"][order(GeneOrder), Gene]
)

zmat <- function(x) {
  x <- as.matrix(x)
  x <- t(scale(t(x)))
  x[!is.finite(x)] <- 0
  x
}

make_matrix <- function(dt, genes) {
  genes <- intersect(genes, names(dt))
  x <- as.matrix(dt[, c("TISI_scaled", genes), with = FALSE])
  m <- zmat(t(x))
  rownames(m) <- c("TISI", genes)
  m
}

heatmap_from_matrix <- function(m, title) {
  Heatmap(
    m,
    name = "Z",
    col = colorRamp2(c(-2, 0, 2), c("#2C7BB6", "white", "#D7191C")),
    column_title = title,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    border = TRUE
  )
}

macro_m <- make_matrix(macro_dt, gene_sets$Macrophages)
fibro_m <- make_matrix(fibro_dt, gene_sets$Fibroblasts)

## Individual reproductions using the original Figure 6G heatmap settings.
pdf(file.path(source_dir, "Figure_6G_reproduced_macrophages.pdf"), width = 14, height = 6, useDingbats = FALSE)
draw(heatmap_from_matrix(macro_m, "Macrophages"), heatmap_legend_side = "right")
dev.off()

pdf(file.path(source_dir, "Figure_6G_reproduced_fibroblasts.pdf"), width = 14, height = 6, useDingbats = FALSE)
draw(heatmap_from_matrix(fibro_m, "Fibroblasts"), heatmap_legend_side = "right")
dev.off()

## Compact two-panel version corresponding to paper Figure 6G.
pdf(file.path(source_dir, "Figure_6G_reproduced_combined.pdf"), width = 14, height = 6, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.42, 5.58), "in"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text("Iron gene expression in tumor neighboring cells", gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(heatmap_from_matrix(macro_m, "Macrophages"), newpage = FALSE, heatmap_legend_side = "right")
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(heatmap_from_matrix(fibro_m, "Fibroblasts"), newpage = FALSE, heatmap_legend_side = "right")
popViewport()
popViewport()
dev.off()

cat("Generated Figure 6G PDFs from source data only.\n")
