suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source_csv <- "Figure_6A_source_data.csv"
output_pdf <- "Figure_6A.pdf"

dt <- fread(source_csv)
required <- c("sample_id", "stage", "tisi_group", "sender", "receiver", "log1p_weight")
stopifnot(all(required %in% names(dt)))

cell_levels <- c(
  "Tumor_cells", "Macrophages", "Endothelial_cells", "Fibroblasts",
  "T_cells", "Dendritic_cells", "B_cells"
)
short_labels <- c(
  Tumor_cells = "Tu", Macrophages = "M", Endothelial_cells = "E",
  Fibroblasts = "F", T_cells = "T", Dendritic_cells = "D", B_cells = "B"
)
row_labels <- c(
  Tumor_cells = "Tumor cells", Macrophages = "Macrophages",
  Endothelial_cells = "Endothelial", Fibroblasts = "Fibroblasts",
  T_cells = "T cells", Dendritic_cells = "DC cells", B_cells = "B cells"
)
sample_order <- c("20221010-1_1", "20221010-2_6", "20220401-2_8", "20221010-2_8")

make_matrix <- function(sid) {
  x <- dt[sample_id == sid]
  m <- matrix(0, length(cell_levels), length(cell_levels),
              dimnames = list(cell_levels, cell_levels))
  m[cbind(match(x$sender, cell_levels), match(x$receiver, cell_levels))] <-
    x$log1p_weight
  m
}

mats <- lapply(sample_order, make_matrix)
global_max <- max(unlist(mats), na.rm = TRUE)
col_fun <- colorRamp2(
  c(0, global_max * 0.15, global_max * 0.5, global_max),
  c("white", "#9EC5E5", "#E07A7A", "#7A0000")
)

stage_text <- c("Primary\nLow", "TAC 7d", "TAC 12d\nHigh", "Recurrent\nLow")
heatmaps <- lapply(seq_along(mats), function(i) {
  m <- mats[[i]]
  colnames(m) <- unname(short_labels[colnames(m)])
  rownames(m) <- unname(row_labels[rownames(m)])
  Heatmap(
    m,
    name = paste0("h", i),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "#4F5D73", lwd = 0.6),
    row_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    show_row_names = i == length(mats),
    show_heatmap_legend = i == 1,
    row_names_side = "right",
    column_title = stage_text[i],
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(title = "log1p(weight)"),
  )
})

pdf(output_pdf, width = 12, height = 3.8, useDingbats = FALSE)
draw(
  Reduce(`+`, heatmaps),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)
grid.text("Receiver", x = 0.39, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"))
grid.text("Sender", x = 0.015, y = 0.50, rot = 90,
          gp = gpar(fontsize = 11, fontface = "bold"))
dev.off()
