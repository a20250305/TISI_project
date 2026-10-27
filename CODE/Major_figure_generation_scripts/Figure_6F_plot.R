suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source_csv <- "Figure_6F_source_data.csv"
manifest_csv <- "Figure_6F__manifest.csv"
output_pdf <- "Figure_6F_full_neighbor_reproduced_from_source_data.pdf"

dt <- fread(source_csv)
manifest <- fread(manifest_csv)
stopifnot(all(c("Gene", "GeneSet") %in% names(dt)))
stopifnot(all(c("GeneSet", "Gene", "GeneOrder") %in% names(manifest)))

gene_set_names <- unique(manifest$GeneSet)
gene_sets <- setNames(
  lapply(gene_set_names, function(x) manifest[GeneSet == x][order(GeneOrder), Gene]),
  gene_set_names
)

celltype_order <- c(
  "Tumor_cells", "Fibroblasts", "Endothelial_cells", "Macrophages",
  "NK_cells", "Dendritic_cells", "T_cells", "B_cells"
)

celltype_labels <- c(
  Tumor_cells = "Tumor cells",
  Fibroblasts = "Fibroblasts",
  Endothelial_cells = "Endothelial cells",
  Macrophages = "Macrophages",
  NK_cells = "NK cells",
  Dendritic_cells = "DC cells",
  T_cells = "T cells",
  B_cells = "B cells"
)

gene_set_labels <- c(
  iron_uptake = "Iron uptake",
  iron_storage = "Iron storage",
  heme_metabolism = "Heme metabolism",
  iron_trafficking = "Iron trafficking",
  iron_TME = "Iron TME"
)

zmat <- function(x) {
  x <- as.matrix(x)
  x <- t(scale(t(x)))
  x[!is.finite(x)] <- 0
  x
}

make_matrix <- function(gene_set) {
  genes <- gene_sets[[gene_set]]
  sub <- dt[match(genes, Gene)]
  sub <- sub[!is.na(Gene)]
  long <- melt(sub, id.vars = c("Gene", "GeneSet"), variable.name = "SampleCellType", value.name = "Expression")
  long[, CellType := sub("^[0-9]{8}-[0-9]+_[0-9]+_", "", SampleCellType)]
  m0 <- dcast(
    long[, .(Expression = mean(Expression, na.rm = TRUE)), by = .(Gene, CellType)],
    Gene ~ CellType,
    value.var = "Expression"
  )
  missing_celltypes <- setdiff(celltype_order, names(m0))
  if (length(missing_celltypes) > 0) {
    for (x in missing_celltypes) m0[, (x) := NA_real_]
  }
  m0 <- m0[, c("Gene", celltype_order), with = FALSE]
  m <- zmat(m0[, -1, with = FALSE])
  rownames(m) <- m0$Gene
  m <- t(m)[celltype_order, , drop = FALSE]
  rownames(m) <- unname(celltype_labels[rownames(m)])
  m
}

make_heatmap <- function(gene_set, show_rows = FALSE) {
  Heatmap(
    make_matrix(gene_set),
    name = "Z",
    col = colorRamp2(c(-2, 0, 2), c("#2C7BB6", "white", "#D7191C")),
    column_title = unname(gene_set_labels[gene_set]),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = show_rows,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 11),
    border = TRUE,
    show_heatmap_legend = FALSE
  )
}

pdf(output_pdf, width = 12, height = 7, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))

gene_set_names <- names(gene_sets)
for (i in seq_along(gene_set_names)) {
  row <- if (i <= 3) 1 else 2
  col <- if (i <= 3) i else i - 3
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  draw(make_heatmap(gene_set_names[i], show_rows = i == 3), newpage = FALSE)
  popViewport()
}

draw(
  Legend(
    title = "Z",
    col_fun = colorRamp2(c(-2, 0, 2), c("#2C7BB6", "white", "#D7191C")),
    at = c(-2, -1, 0, 1, 2),
    direction = "vertical"
  ),
  x = unit(0.91, "npc"),
  y = unit(0.28, "npc")
)

popViewport()
dev.off()

cat("Generated:", output_pdf, "\n")
