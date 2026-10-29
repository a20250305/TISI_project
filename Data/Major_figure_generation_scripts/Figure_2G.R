suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_2G_source_data.csv"
output_pdf <- "Figure_2G_reproduced.pdf"

d <- fread(source_csv)
stopifnot(all(c("CancerType", "Dataset", "n") %in% names(d)))

d[, n := as.numeric(n)]
d[, FigureLabel := fifelse(CancerType == "EADC", "EAC", CancerType)]

stopifnot(
  sum(d$n) == 1380,
  uniqueN(d$Dataset) == 15,
  uniqueN(d$CancerType) == 4
)

plot_order <- c("BRCA", "EAC", "OV", "READ")
ct <- d[, .(n = sum(n)), by = FigureLabel]
ct[, FigureLabel := factor(FigureLabel, levels = plot_order)]
setorder(ct, FigureLabel)
ct[, ymax := cumsum(n)]
ct[, ymin := shift(ymax, fill = 0)]

palette <- c(
  BRCA = "#9F85AF",
  READ = "#B9E37A",
  OV   = "#67BFAB",
  EAC  = "#399BD1"
)
legend_order <- c("BRCA", "READ", "OV", "EAC")

p <- ggplot() +
  geom_rect(
    data = ct,
    aes(ymin = ymin, ymax = ymax, xmin = 2.2, xmax = 3, fill = FigureLabel),
    color = "white",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0,
    y = sum(ct$n) / 2,
    label = "1,380 tumor samples\n15 datasets\nPre- & post-chemo",
    size = 4.8,
    lineheight = 0.95
  ) +
  coord_polar(theta = "y", clip = "off") +
  xlim(0, 3) +
  scale_fill_manual(values = palette, breaks = legend_order) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "right",
    plot.margin = margin(4, 12, 4, 4)
  )

ggsave(
  output_pdf,
  p,
  width = 5.8,
  height = 4.8,
  units = "in",
  device = "pdf"
)
