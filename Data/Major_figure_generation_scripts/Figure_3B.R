suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3B_source_data.csv"
output_pdf <- "Figure_3B_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("perturbation", "n_cells", "deltaMean", "direction") %in% names(dt)))
stopifnot(all(is.finite(dt$n_cells)), all(dt$n_cells > 0))
stopifnot(all(is.finite(dt$deltaMean)))

p <- ggplot(dt, aes(n_cells, deltaMean)) +
  geom_point(aes(color = direction), alpha = 0.7, size = 1.2) +
  scale_x_log10() +
  scale_color_manual(values = c(
    "TISI increase" = "#B2182B",
    "TISI decrease" = "#2166AC",
    "Neutral" = "grey70"
  )) +
  theme_classic(base_size = 14) +
  labs(
    title = "TISI effect size vs number of cells",
    x = "Number of cells (log10)",
    y = "Δ TISI (scaled)",
    color = NULL
  )

ggsave(output_pdf, p, width = 7, height = 6, units = "in", device = "pdf")
