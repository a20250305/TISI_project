suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3D_source_data.csv"
output_pdf <- "Figure_3D_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c(
  "Sample", "Treatment", "Genotype", "Group",
  "TISI_raw", "TISI_scaled", "order_index"
) %in% names(dt)))
stopifnot(nrow(dt) == 12, all(is.finite(dt$TISI_scaled)))

group_order <- c(
  "WT Vehicle", "KO Vehicle",
  "WT Gemcitabine", "KO Gemcitabine"
)
dt$Group <- factor(dt$Group, levels = group_order)
dt$x <- match(dt$Group, group_order)

p_vehicle <- t.test(
  TISI_scaled ~ Genotype,
  data = dt[dt$Treatment == "Vehicle", ],
  var.equal = TRUE
)$p.value
p_gemcitabine <- t.test(
  TISI_scaled ~ Genotype,
  data = dt[dt$Treatment == "Gemcitabine", ],
  var.equal = TRUE
)$p.value

star <- function(p) {
  if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

fill_colors <- c(
  "WT Vehicle" = "#2C7BB6",
  "KO Vehicle" = "#74ADD1",
  "WT Gemcitabine" = "#FDAE61",
  "KO Gemcitabine" = "#D7191C"
)

p <- ggplot(dt, aes(x = x, y = TISI_scaled, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    colour = "#333333",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.06,
    height = 0,
    size = 2.2,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("WT", "KO", "WT", "KO"),
    limits = c(0.45, 4.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.3, 0.6, 0.9, 1.2),
    limits = c(-0.18, 1.2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = fill_colors) +
  annotate("segment", x = 0.75, xend = 2.25, y = -0.08, yend = -0.08) +
  annotate("segment", x = 2.75, xend = 4.25, y = -0.08, yend = -0.08) +
  annotate("text", x = 1.5, y = -0.15, label = "Vehicle", size = 5) +
  annotate("text", x = 3.5, y = -0.15, label = "Gemcitabine", size = 5) +
  annotate("segment", x = 1, xend = 1, y = 1.04, yend = 1.07) +
  annotate("segment", x = 1, xend = 2, y = 1.07, yend = 1.07) +
  annotate("segment", x = 2, xend = 2, y = 1.04, yend = 1.07) +
  annotate("text", x = 1.5, y = 1.10, label = star(p_vehicle), size = 5) +
  annotate("segment", x = 3, xend = 3, y = 1.04, yend = 1.07) +
  annotate("segment", x = 3, xend = 4, y = 1.07, yend = 1.07) +
  annotate("segment", x = 4, xend = 4, y = 1.04, yend = 1.07) +
  annotate("text", x = 3.5, y = 1.10, label = star(p_gemcitabine), size = 5) +
  labs(x = NULL, y = "TISI") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 15),
    plot.margin = margin(8, 8, 38, 8)
  )

ggsave(output_pdf, p, width = 3.7, height = 3.6, units = "in", device = "pdf")
