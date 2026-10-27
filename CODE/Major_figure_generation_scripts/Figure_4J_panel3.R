suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_4J_panel3_source_data.csv"
output_pdf <- "Figure_4J_panel3_TISI_LCN2.pdf"

dt <- fread(source_csv)
stopifnot(all(c("sample_id", "group", "TISI_scaled", "LCN2") %in% names(dt)))

r_value <- cor(dt$TISI_scaled, dt$LCN2, method = "pearson")
p_value <- cor.test(
  dt$TISI_scaled,
  dt$LCN2,
  method = "spearman"
)$p.value
label <- sprintf("r=%.2f\nP=%.2g", r_value, p_value)

p <- ggplot(dt, aes(TISI_scaled, LCN2)) +
  geom_point(color = "#2C7BB6", size = 2) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = label,
    hjust = -0.1,
    vjust = 1.2,
    size = 3.5
  ) +
  labs(x = "TISI", y = "LCN2") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.4
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(
  output_pdf,
  p,
  width = 2.8,
  height = 2.2,
  units = "in",
  device = "pdf"
)
