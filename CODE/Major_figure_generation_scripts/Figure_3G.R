suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3G_source_data.csv"
output_pdf <- "Figure_3G_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE)
stopifnot(all(c("drug", "rho", "FDR") %in% names(dt)))
stopifnot(nrow(dt) == 14, all(is.finite(dt$rho)), all(is.finite(dt$FDR)))

dt$neglogFDR <- -log10(dt$FDR)
dt$plot_drug <- dt$drug
dt$plot_drug[dt$plot_drug == "Bleomycin (50 uM)"] <- "Bleomycin"
dt$plot_drug <- factor(
  dt$plot_drug,
  levels = dt$plot_drug[order(dt$rho)]
)

p <- ggplot(dt, aes(x = rho, y = plot_drug)) +
  geom_segment(
    aes(x = 0, xend = rho, yend = plot_drug),
    colour = "grey60",
    linewidth = 0.6
  ) +
  geom_point(
    aes(size = neglogFDR, colour = rho),
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = c(-0.1, 0, 0.1, 0.2),
    limits = c(-0.16, 0.24),
    expand = c(0, 0)
  ) +
  scale_colour_gradient2(
    low = "#4C72B0",
    mid = "white",
    high = "#9E1B15",
    midpoint = 0,
    limits = c(-0.16, 0.24),
    breaks = c(-0.1, 0, 0.1, 0.2),
    name = "rho"
  ) +
  scale_size_continuous(
    breaks = c(2.5, 5, 7.5),
    range = c(2, 7),
    name = "-log10(FDR)"
  ) +
  guides(
    colour = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    x = "Spearman correlation (TISI vs IC50)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )

ggsave(
  output_pdf,
  p,
  width = 5.8,
  height = 3.8,
  units = "in",
  device = "pdf"
)
