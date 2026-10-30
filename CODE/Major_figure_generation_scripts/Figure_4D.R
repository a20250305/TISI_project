suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_4D_source_data.csv"
output_pdf <- "Figure_4D.pdf"

dt <- fread(source_csv)
required_columns <- c("cell_id", "sample", "stage", "TISI")
stopifnot(all(required_columns %in% names(dt)))

dt[, TISI := as.numeric(TISI)]
dt[, stage := factor(
  stage,
  levels = c("Naive", "MRD", "Relapse")
)]
dt <- dt[is.finite(TISI) & !is.na(stage)]

stage_colors <- c(
  Naive = "#2C7FB8",
  MRD = "#D7191C",
  Relapse = "#4DAF4A"
)

p <- ggplot(dt, aes(x = stage, y = TISI, fill = stage)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    color = "#2B2B2B",
    linewidth = 0.55
  ) +
  scale_fill_manual(values = stage_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "TISI score") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12
    ),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  filename = output_pdf,
  plot = p,
  width = 3.4,
  height = 3.4,
  units = "in",
  device = "pdf"
)
