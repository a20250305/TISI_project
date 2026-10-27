suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_4E_source_data.csv"
output_pdf <- "PDX98B_Hallmark_bubble_NatureStyle.pdf"

dt <- fread(source_csv)
required_columns <- c("Signature", "GeneRatio", "Count", "p_adjust", "plot_order")
stopifnot(all(required_columns %in% names(dt)))

dt[, GeneRatio := as.numeric(GeneRatio)]
dt[, Count := as.numeric(Count)]
dt[, p_adjust := as.numeric(p_adjust)]
dt[, plot_order := as.numeric(plot_order)]
dt[, logP := -log10(p_adjust)]
dt <- dt[is.finite(GeneRatio) & is.finite(Count) & is.finite(logP)]
setorder(dt, plot_order)
dt[, Signature := factor(Signature, levels = rev(Signature))]

p <- ggplot(
  dt,
  aes(
    x = GeneRatio,
    y = Signature,
    size = Count,
    color = logP
  )
) +
  geom_point(alpha = 0.9, stroke = 0.3) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B") +
  scale_size(range = c(3, 9)) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    color = "-log10(adj P)",
    size = "Gene count"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = output_pdf,
  plot = p,
  width = 6.5,
  height = 6,
  units = "in",
  device = "pdf"
)
