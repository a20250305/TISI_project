suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_5J_source_data.csv"
output_pdf <- "Figure_5J.pdf"

dt <- fread(source_csv)
stopifnot(all(c(
  "sample_id", "group", "x_label", "spot_id", "TISI_scaled"
) %in% names(dt)))

dt[, x_label := factor(
  x_label,
  levels = c("P17-Pre", "P17-Post", "P20-Pre", "P20-Post")
)]
dt[, group := factor(group, levels = c("pre", "post"))]

p <- ggplot(dt, aes(x_label, TISI_scaled, fill = group)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c(pre = "#1F78B4", post = "#E31A1C")) +
  scale_y_continuous(
    breaks = c(0, 0.4, 0.8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(x = NULL, y = "TISI") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12)
  )

ggsave(output_pdf, p, width = 3.8, height = 2.7, units = "in", device = "pdf")
