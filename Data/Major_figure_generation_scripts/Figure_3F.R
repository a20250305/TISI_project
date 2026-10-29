suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source_csv <- "Figure_3F_source.csv"
output_pdf <- "Figure_3F_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE)
stopifnot(all(c(
  "cell_line", "drug", "genotype", "replicate", "viability_percent"
) %in% names(dt)))

drug_order <- c("Gem", "Eto", "Iri", "Pac")
genotype_order <- c("WT", "MED12 KO")
dt$drug <- factor(dt$drug, levels = drug_order)
dt$genotype <- factor(dt$genotype, levels = genotype_order)
dt$x <- as.numeric(dt$drug) + ifelse(dt$genotype == "WT", -0.17, 0.17)

summary_dt <- aggregate(
  viability_percent ~ cell_line + drug + genotype + x,
  data = dt,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)
summary_dt$mean <- summary_dt$viability_percent[, "mean"]
summary_dt$sd <- summary_dt$viability_percent[, "sd"]
summary_dt$viability_percent <- NULL

stars <- function(p) ifelse(p < 0.001, "**", ifelse(p < 0.05, "*", "ns"))

make_panel <- function(cell, y_breaks, y_limit, show_y = FALSE) {
  d <- dt[dt$cell_line == cell, ]
  s <- summary_dt[summary_dt$cell_line == cell, ]

  p_values <- sapply(drug_order, function(drug) {
    x <- d[d$drug == drug, ]
    t.test(viability_percent ~ genotype, data = x)$p.value
  })

  bracket <- data.frame(
    drug = drug_order,
    x1 = seq_along(drug_order) - 0.17,
    x2 = seq_along(drug_order) + 0.17,
    y = sapply(drug_order, function(drug) {
      max(s$mean[s$drug == drug] + s$sd[s$drug == drug]) * 1.15
    }),
    label = stars(p_values)
  )

  p <- ggplot() +
    geom_col(
      data = s,
      aes(x = x, y = mean, color = genotype),
      fill = "white",
      width = 0.30,
      linewidth = 0.7
    ) +
    geom_errorbar(
      data = s,
      aes(x = x, ymin = mean - sd, ymax = mean + sd, color = genotype),
      width = 0.08,
      linewidth = 0.7
    ) +
    geom_point(
      data = d,
      aes(x = x, y = viability_percent, color = genotype),
      size = 1.8,
      position = position_jitter(width = 0.025, height = 0)
    ) +
    geom_segment(
      data = bracket,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket,
      aes(x = x1, xend = x1, y = y - diff(y_limit) * 0.018, yend = y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = bracket,
      aes(x = x2, xend = x2, y = y - diff(y_limit) * 0.018, yend = y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_text(
      data = bracket,
      aes(x = (x1 + x2) / 2, y = y + diff(y_limit) * 0.018, label = label),
      inherit.aes = FALSE,
      size = 4.2,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = seq_along(drug_order),
      labels = c("Gem.", "Eto.", "Iri.", "Pac."),
      limits = c(0.45, 4.55),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_limit,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = c("WT" = "black", "MED12 KO" = "#E31A1C")) +
    labs(
      title = cell,
      x = NULL,
      y = if (show_y) "Cell viability (% untreated)" else NULL,
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5)
    )

  p
}

p1 <- make_panel("Mia PaCa-2", c(0, 20, 40, 60), c(0, 70), TRUE)
p2 <- make_panel("MDA-MB-468", c(0, 20, 40, 60, 80, 100), c(0, 120))
p3 <- make_panel("MDA-MB-231", c(0, 20, 40, 60, 80), c(0, 95))

figure <- (p1 + p2 + p3) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(output_pdf, figure, width = 8.2, height = 3.7, units = "in", device = "pdf")
