suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source_csv <- "Figure_5E_source_data.csv"
output_pdf <- "Figure_5E.pdf"

dt <- fread(source_csv)
stopifnot(all(c("Niche", "Feature", "rho") %in% names(dt)))
dt[, Niche := factor(Niche, levels = paste0("Niche_", 1:13))]
dt[, Feature := as.character(Feature)]

rho <- dcast(dt, Feature ~ Niche, value.var = "rho")
mat <- as.matrix(rho[, -1])
rownames(mat) <- rho$Feature
storage.mode(mat) <- "numeric"

feature_hclust <- hclust(dist(mat))
rho_col_fun <- colorRamp2(
  c(-0.5, 0, 0.5),
  c("#2166AC", "white", "#B2182B")
)

pdf(output_pdf, width = 10, height = 8)
heatmap <- Heatmap(
  t(mat),
  name = "rho",
  col = rho_col_fun,
  cluster_rows = TRUE,
  cluster_columns = feature_hclust,
  row_dend_width = unit(20, "mm"),
  column_dend_height = unit(22, "mm"),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "rho")
)
draw(heatmap)
dev.off()
