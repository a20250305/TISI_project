suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3B_source_data.csv"
output_pdf <- "ALL_deltaTISI_vs_nCells_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("perturbation", "n_cells", "deltaMean", "direction") %in% names(dt)))
stopifnot(all(is.finite(dt$n_cells)), all(dt$n_cells > 0), all(is.finite(dt$deltaMean)))

p_all_cells <- ggplot(dt, aes(x = n_cells, y = deltaMean)) +
  geom_point(
    aes(color = direction),
    alpha = 0.7,
    size = 1.2
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "TISI increase" = "#B2182B",
      "TISI decrease" = "#2166AC",
      "Neutral"       = "grey70"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "TISI effect size vs number of cells",
    x = "Number of cells (log10)",
    y = expression(Delta~TISI~"(scaled)"),
    color = NULL
  )

ggsave(output_pdf, p_all_cells, width = 7, height = 6, units = "in")
