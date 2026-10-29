suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

violin_data <- fread("Figure_3I_violin_source_data.csv")
response_data <- fread("Figure_3I_response_source_data.csv")

stopifnot(
  all(c("PDX", "sample", "TISI") %in% names(violin_data)),
  all(c("PDX", "baseline_sample", "TISI", "tumor_size_change_pct") %in% names(response_data))
)

pdx_order <- c("PDX27", "PDX235", "PDX626")
violin_data[, PDX := factor(PDX, levels = pdx_order)]
response_data[, PDX := factor(PDX, levels = pdx_order)]

violin_colors <- c(
  PDX27 = "#2166AC",
  PDX235 = "#B2182B",
  PDX626 = "#B2182B"
)

p_left <- ggplot(violin_data, aes(PDX, TISI, fill = PDX)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "#333333") +
  scale_fill_manual(values = violin_colors) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "TISI") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.margin = margin(3, 3, 3, 3)
  )

response_data[, label := ifelse(
  tumor_size_change_pct >= 0,
  sprintf("+%g%%", tumor_size_change_pct),
  sprintf("%g%%", tumor_size_change_pct)
)]
response_data[, label_hjust := ifelse(tumor_size_change_pct >= 0, 1.1, -0.1)]

p_right <- ggplot(
  response_data,
  aes(x = PDX, y = tumor_size_change_pct, fill = TISI)
) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey40") +
  geom_text(
    aes(y = 0, label = label, hjust = label_hjust),
    size = 3.2
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-110, 125),
    breaks = c(-100, -50, 0, 50, 100),
    expand = c(0, 0)
  ) +
  scale_fill_gradient(
    low = "#2166AC",
    high = "#B2182B",
    limits = range(response_data$TISI),
    breaks = range(response_data$TISI),
    labels = c("low", "high"),
    name = "TISI",
    guide = guide_colorbar(
      ticks = FALSE,
      frame.colour = "black",
      barheight = grid::unit(0.85, "cm"),
      barwidth = grid::unit(0.16, "cm")
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "OV with high TISI has poor response"
  ) +
  theme_classic(base_size = 9) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.5
    ),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )

figure <- p_left + p_right + plot_layout(widths = c(0.82, 1.90))

ggsave(
  "Figure_3I_reproduced.pdf",
  figure,
  width = 5.8,
  height = 2.35,
  units = "in",
  device = "pdf"
)
