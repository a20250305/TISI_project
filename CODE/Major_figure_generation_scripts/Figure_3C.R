suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3C_source_data.csv"
output_pdf <- "CELL_TISI_boxplot_distinct_muted_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("cell_barcode", "perturbation", "TISI", "group", "order_index") %in% names(dt)))
stopifnot(all(is.finite(dt$TISI)))

order_table <- unique(dt[c("perturbation", "order_index")])
order_table <- order_table[order(order_table$order_index), , drop = FALSE]

if (nrow(order_table) == 0L ||
    !identical(order_table$order_index, seq_len(nrow(order_table)))) {
  stop("order_index must be a complete sequence starting at 1.")
}

if (anyDuplicated(order_table$perturbation) > 0L) {
  stop("Each perturbation must have one unique order_index.")
}

order_levels <- order_table$perturbation
dt$perturbation <- factor(dt$perturbation, levels = order_levels)

palette_okabe <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#999999", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#8DA0CB",
  "#FC8D62", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
  "#4D4D4D"
)

all_perts <- levels(dt$perturbation)
cols_B <- setNames(
  rep(palette_okabe, length.out = length(all_perts)),
  all_perts
)

if ("non-targeting" %in% names(cols_B)) {
  cols_B["non-targeting"] <- "grey60"
}

p_box_B <- ggplot(
  dt,
  aes(x = perturbation, y = TISI, fill = perturbation)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  scale_fill_manual(values = cols_B) +
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

ggsave(output_pdf, p_box_B, width = 7, height = 5, units = "in")
