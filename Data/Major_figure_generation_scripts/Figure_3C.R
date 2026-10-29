suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3C_source_data.csv"
output_pdf <- "Figure_3C_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("cell_barcode", "perturbation", "TISI", "group", "order_index") %in% names(dt)))
stopifnot(all(is.finite(dt$TISI)))

order_table <- unique(dt[c("perturbation", "order_index")])
order_table <- order_table[order(order_table$order_index), , drop = FALSE]
stopifnot(identical(order_table$order_index, seq_len(nrow(order_table))))
stopifnot(!anyDuplicated(order_table$perturbation))

dt$perturbation <- factor(
  dt$perturbation,
  levels = order_table$perturbation
)

palette_okabe <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#999999", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#8DA0CB",
  "#FC8D62", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
  "#4D4D4D"
)

plot_colors <- setNames(
  rep(palette_okabe, length.out = nlevels(dt$perturbation)),
  levels(dt$perturbation)
)
if ("non-targeting" %in% names(plot_colors)) {
  plot_colors["non-targeting"] <- "grey60"
}

p <- ggplot(dt, aes(perturbation, TISI, fill = perturbation)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = plot_colors) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Cell-level TISI (distinct colors per perturbation)",
    subtitle = "Muted palette (no rainbow), box=IQR; whiskers=1.5×IQR",
    x = NULL,
    y = "TISI score"
  )

ggsave(output_pdf, p, width = 7, height = 5, units = "in", device = "pdf")
