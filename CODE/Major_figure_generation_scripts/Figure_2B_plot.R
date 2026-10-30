suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_2B_source_data.csv"
output_pdf <- "Figure_2B_reproduced.pdf"

dt <- fread(source_csv)
stopifnot(all(c("Dataset", "Method", "AUC") %in% names(dt)))
dt[, AUC := as.numeric(AUC)]
stopifnot(length(unique(dt$Method)) == 32L)

dataset_order <- c(
  "GSE146590", "GSE234387", "GSE276238", "GSE280381_1A",
  "GSE280381_1B", "GSE280381_2A", "GSE280381_2B", "GSE293456",
  "GSE297393"
)
dt[, Dataset := factor(Dataset, levels = dataset_order)]

method_order <- dt[!is.na(AUC), .(mean_AUC = mean(AUC)), by = Method][
  order(mean_AUC), Method
]
method_order <- c(setdiff(method_order, "TISI"), "TISI")
dt[, Method := factor(Method, levels = method_order)]

method_labels <- setNames(method_order, method_order)
method_labels[c("Reactome_CS", "Reactome_OSIS", "CSS20")] <-
  c("REACT_CS", "REAC_OSIS", "CSS")
display_levels <- unname(method_labels[method_order])
dt[, Method_label := factor(
  method_labels[as.character(Method)], levels = display_levels
)]

bar <- dt[!is.na(AUC), .(AUC = mean(AUC)), by = Method]
bar[, method_name := as.character(Method)]
bar[, Method_label := factor(
  method_labels[method_name], levels = display_levels
)]

p <- ggplot() +
  geom_col(
    data = bar[method_name != "TISI"],
    aes(x = Method_label, y = AUC),
    fill = "white", color = "grey55", linewidth = 0.7, width = 0.62
  ) +
  geom_col(
    data = bar[method_name == "TISI"],
    aes(x = Method_label, y = AUC),
    fill = "#F6B6B6", color = "#E41A1C", linewidth = 0.8, width = 0.62
  ) +
  geom_point(
    data = dt[!is.na(AUC)],
    aes(x = Method_label, y = AUC, color = Dataset),
    size = 1.8, alpha = 0.95,
    position = position_jitter(width = 0.3, height = 0.02, seed = 1)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_color_viridis_d(option = "C", end = 0.9, drop = FALSE) +
  labs(title = "Single cell RNA-seq", x = NULL, y = "AUC", color = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 22, 30)
  )

ggsave(output_pdf, p, width = 8.5, height = 4.2, units = "in", device = "pdf")
