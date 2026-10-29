suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_2E_source_data.csv"
output_pdf <- "Figure_2E_reproduced.pdf"

d <- fread(source_csv)
stopifnot(all(c("group", "replicate", "value") %in% names(d)))
d[, group := factor(group, levels = c("Control", "TIS"))]
d[, value := as.numeric(value)]

summary_d <- d[, .(mean = mean(value), sd = sd(value)), by = group]
p_value <- t.test(value ~ group, data = d)$p.value
label <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else if (p_value < 0.05) "*" else "ns"

p <- ggplot(summary_d, aes(group, mean)) +
  geom_col(width = 0.55, fill = "white", color = "black", linewidth = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
    width = 0.10, linewidth = 0.7
  ) +
  geom_jitter(
    data = d,
    aes(group, value, color = group),
    width = 0.07, height = 0, size = 2.3,
    show.legend = FALSE
  ) +
  annotate("text", x = 1.5, y = 3820, label = label, size = 6) +
  scale_color_manual(values = c(Control = "#377EB8", TIS = "#E15759")) +
  scale_y_continuous(
    breaks = c(0, 1000, 2000, 3000),
    labels = c("0", "1K", "2K", "3K"),
    limits = c(0, 4000),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "MFI") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(output_pdf, p, width = 2.6, height = 3.2, units = "in", device = "pdf")
