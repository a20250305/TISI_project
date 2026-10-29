suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_2F_source_data.csv"
output_pdf <- "Figure_2F_reproduced.pdf"

d <- fread(source_csv)
stopifnot(all(c("Dataset", "cellname", "label", "TISI") %in% names(d)))

d[, TISI := as.numeric(TISI)]
d[, group := factor(
  label,
  levels = c("non-Sene", "Sene"),
  labels = c("Control", "TIS")
)]

p_value <- t.test(TISI ~ group, data = d)$p.value

p <- ggplot(d, aes(group, TISI, fill = group)) +
  geom_violin(
    trim = FALSE,
    width = 0.9,
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    width = 0.20,
    fill = "white",
    color = "black",
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  annotate("text", x = 1.5, y = 1.02, label = "**", size = 5.5) +
  scale_fill_manual(values = c(Control = "#2C7BB6", TIS = "#D7191C")) +
  scale_y_continuous(
    breaks = c(0, 0.3, 0.6, 0.9),
    limits = c(-0.1, 1.1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "TISI") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.6
    ),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  output_pdf,
  p,
  width = 2.6,
  height = 3.2,
  units = "in",
  device = "pdf"
)

message("Welch t-test P = ", format.pval(p_value, digits = 4))
