suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

source_csv <- "Figure_6B_source_data.csv"
output_pdf <- "Figure_6B.pdf"

dt <- fread(source_csv)
required <- c("sender", "receiver", "lr", "prob", "selected_for_panel_B")
stopifnot(all(required %in% names(dt)))

mat <- dcast(
  dt[selected_for_panel_B == TRUE],
  sender ~ receiver,
  value.var = "prob",
  fun.aggregate = sum,
  fill = 0
)
rn <- mat$sender
mat$sender <- NULL
mat <- as.matrix(mat)
rownames(mat) <- rn

active <- rowSums(mat) + colSums(mat) > 0
mat <- mat[active, active, drop = FALSE]

cell_labels <- c(
  Tumor_cells = "Tumor cells", Macrophages = "Macrophages",
  Fibroblasts = "Fibroblast", B_cells = "B cells",
  T_cells = "T cells", Dendritic_cells = "DC cells"
)
cell_cols <- c(
  Tumor_cells = "#E64B35",
  Macrophages = "#4DBBD5",
  Fibroblasts = "#00A087",
  B_cells = "#3C5488",
  T_cells = "#F39B7F",
  Dendritic_cells = "#8491B4"
)
cell_cols <- cell_cols[rownames(mat)]
names(cell_cols) <- unname(cell_labels[rownames(mat)])
dimnames(mat) <- list(unname(cell_labels[rownames(mat)]),
                      unname(cell_labels[colnames(mat)]))

pdf(output_pdf, width = 7.4, height = 7.4, useDingbats = FALSE)
par(mar = c(5, 1, 3, 1))
circos.clear()
circos.par(
  track.margin = c(0.01, 0.01),
  start.degree = 90,
  gap.degree = 5,
  track.height = 0.12
)
chordDiagram(
  mat,
  directional = 1,
  direction.type = "arrows",
  transparency = 0.25,
  grid.col = cell_cols,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    circos.text(
      mean(get.cell.meta.data("xlim")),
      get.cell.meta.data("ylim")[1] + 0.1,
      sector,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.8,
      font = 2
    )
  },
  bg.border = NA
)
title("Iron-related ligand-receptor pairs", cex.main = 1.15, font.main = 2)
mtext("TAC 12d", side = 1, line = 0.5, cex = 1.15, font = 2)
dev.off()
circos.clear()
