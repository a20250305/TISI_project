suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source_csv <- "Figure_5D_source_data.csv"
output_pdf <- "Figure_5D.pdf"

dt <- fread(source_csv)
stopifnot(all(c("Niche", "Major", "Fraction") %in% names(dt)))
dt[, Fraction := as.numeric(Fraction)]
dt <- dt[is.finite(Fraction)]

major_order <- c(
  "B_cells", "Dendritic_cells", "Endothelial_cells", "Fibroblasts",
  "Macrophages", "NK_cells", "T_cells", "Tumor_cells"
)

palette <- c(
  B_cells = "#F4D03F",
  Dendritic_cells = "#7F8C8D",
  Endothelial_cells = "#16A085",
  Fibroblasts = "#C0392B",
  Macrophages = "#D68910",
  NK_cells = "#5DADE2",
  T_cells = "#2471A3",
  Tumor_cells = "#8E44AD"
)

dt[, Major := factor(Major, levels = major_order)]
wide <- dcast(dt, Niche ~ Major, value.var = "Fraction", fill = 0)
wide <- wide[, c("Niche", major_order), with = FALSE]
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$Niche

hc <- hclust(dist(mat), method = "ward.D2")
hc$height <- log1p(hc$height)
niche_order <- rownames(mat)[hc$order]
dt[, Niche := factor(Niche, levels = niche_order)]

get_xy <- function(node) {
  if (node < 0) {
    return(c(which(niche_order == rownames(mat)[-node]), 0))
  }
  left <- get_xy(hc$merge[node, 1])
  right <- get_xy(hc$merge[node, 2])
  c(mean(c(left[1], right[1])), hc$height[node])
}

segments <- rbindlist(lapply(seq_len(nrow(hc$merge)), function(i) {
  left <- get_xy(hc$merge[i, 1])
  right <- get_xy(hc$merge[i, 2])
  rbind(
    data.table(x = left[1], y = left[2], xend = left[1], yend = hc$height[i]),
    data.table(x = right[1], y = right[2], xend = right[1], yend = hc$height[i]),
    data.table(x = left[1], y = hc$height[i], xend = right[1], yend = hc$height[i])
  )
}))

p_dend <- ggplot(segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4) +
  scale_x_continuous(limits = c(1, length(niche_order)), expand = c(0, 0)) +
  theme_void()

p_bar <- ggplot(dt, aes(Niche, Fraction, fill = Major)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  xlab("Spatial niche") +
  ylab("Cell composition")

ggsave(
  output_pdf,
  p_dend / p_bar + plot_layout(heights = c(1, 3)),
  width = 10,
  height = 6,
  units = "in",
  device = "pdf"
)
