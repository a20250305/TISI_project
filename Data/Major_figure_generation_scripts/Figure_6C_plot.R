suppressPackageStartupMessages({
  library(data.table)
  library(CellChat)
})

source_csv <- "Figure_6C_source_data.csv"
output_pdf <- "Figure_6C.pdf"

dt <- fread(source_csv)
required <- c("sender", "receiver", "lr", "prob", "panel_order")
stopifnot(all(required %in% names(dt)))

cell_cols <- c(
  Tumor_cells = "#E64B35",
  Macrophages = "#4DBBD5",
  Fibroblasts = "#00A087",
  B_cells = "#3C5488",
  T_cells = "#F39B7F",
  Dendritic_cells = "#8491B4"
)
cell_labels <- c(
  Tumor_cells = "Tumor cells", Macrophages = "Macrophages",
  Fibroblasts = "Fibroblasts", B_cells = "B cells",
  T_cells = "T cells", Dendritic_cells = "DC cells"
)
lr_order <- c("LCN2_SLC22A17", "FTH1_TFRC", "LTF_LRP1", "HFE_TFRC")
lr_titles <- c(
  LCN2_SLC22A17 = "LCN2-SLC22A17",
  FTH1_TFRC = "FTH1-TFRC",
  LTF_LRP1 = "LTF-LRP1",
  HFE_TFRC = "HFE-TFRC"
)

make_matrix <- function(lr_name) {
  x <- dt[lr == lr_name]
  mat <- dcast(x, sender ~ receiver, value.var = "prob", fill = 0)
  rn <- mat$sender
  mat$sender <- NULL
  mat <- as.matrix(mat)
  rownames(mat) <- rn
  all_cells <- union(rownames(mat), colnames(mat))
  full <- matrix(0, length(all_cells), length(all_cells),
                 dimnames = list(all_cells, all_cells))
  full[rownames(mat), colnames(mat)] <- mat
  active <- rowSums(full) + colSums(full) > 0
  full <- full[active, active, drop = FALSE]
  dimnames(full) <- list(unname(cell_labels[rownames(full)]),
                         unname(cell_labels[colnames(full)]))
  full
}

pdf(output_pdf, width = 7.2, height = 7.2, useDingbats = FALSE)
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
for (lr_name in lr_order) {
  mat <- make_matrix(lr_name)
  netVisual_circle(
    mat,
    vertex.weight = rowSums(mat),
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = setNames(cell_cols[names(cell_labels)[match(rownames(mat), cell_labels)]],
                         rownames(mat)),
    edge.width.max = 28,
    vertex.label.cex = 1.1
  )
  title(main = lr_titles[[lr_name]], cex.main = 1.15, font.main = 2)
}
dev.off()
