suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

source_csv <- "Figure_5K_source_data.csv"
output_pdf <- "Figure_5K.pdf"

dt <- fread(source_csv)
stopifnot(all(c(
  "sample_id", "feature", "rho", "p", "neglog10_p", "included"
) %in% names(dt)))
dt <- dt[included == TRUE]

sample_levels <- c("P14A", "P16A", "P17A", "P18A", "P19A", "P20A")
feature_levels <- c(
  "B_cell", "DC", "Endothelium", "Epithelium", "Fibroblast",
  "Lymphocyte", "Macrophage", "Neutrophil", "Plasma"
)

dt[, sample_label := factor(
  sub("A$", "", sample_id),
  levels = sub("A$", "", sample_levels)
)]
dt[, feature := factor(feature, levels = feature_levels)]

p <- ggplot(dt, aes(feature, sample_label)) +
  geom_point(
    aes(size = neglog10_p, fill = rho),
    shape = 21, color = "grey30", stroke = 0.35
  ) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    limits = c(-0.4, 0.4), oob = squish, name = "Neighbor rho"
  ) +
  scale_size(
    range = c(1.5, 6), limits = c(0, 55),
    breaks = c(10, 20, 30, 40, 50), name = "-log10 P"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.spacing.y = grid::unit(0.08, "cm")
  )

ggsave(output_pdf, p, width = 5.2, height = 3.2, units = "in", device = "pdf")
