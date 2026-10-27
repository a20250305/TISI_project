suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_4L_panel1_source_data.csv"
output_pdf <- "Figure_4L_panel1_TISI.pdf"

dt <- fread(source_csv)
stopifnot(all(c("sample_id", "group", "score", "included") %in% names(dt)))
dt <- dt[included == TRUE]
dt[, group := factor(group, levels = c("pre", "post"))]

p_value <- with(dt, wilcox.test(
  score[group == "post"],
  score[group == "pre"]
)$p.value)
label <- ifelse(
  p_value < 0.001, "***",
  ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns"))
)

set.seed(260305)
p <- ggplot(dt, aes(group, score, fill = group)) +
  geom_boxplot(
    width = 0.6,
    color = "black",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(width = 0.13),
    color = "black",
    size = 1.1
  ) +
  annotate("text", x = 1.5, y = Inf, label = label,
           vjust = 1.2, size = 4.5) +
  scale_fill_manual(values = c(pre = "#2C7BB6", post = "#D7191C")) +
  labs(x = NULL, y = "TISI") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.35
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(output_pdf, p, width = 2.2, height = 2.2, units = "in", device = "pdf")
