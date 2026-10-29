suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_5C_source_data.csv"
output_pdf <- "Figure_5C.pdf"

sample_order <- c(
  "20220401-1_3", "20220401-2_8", "20221010-2_7", "20221010-3_11"
)
sample_labels <- c(
  "20220401-1_3" = "MRD-1",
  "20220401-2_8" = "MRD-2",
  "20221010-2_7" = "MRD-3",
  "20221010-3_11" = "MRD-4"
)
niche_order <- paste0("Niche_", 1:13)

dt <- fread(source_csv)
self_dt <- dt[grepl("^Self_[0-9]+$", Variable)]
self_dt[, niche_order := as.integer(sub("^Self_", "", Variable))]
setorder(self_dt, niche_order)
stopifnot(nrow(self_dt) == 13L)

value_order <- rev(sample_order)
plot_dt <- rbindlist(lapply(seq_along(sample_order), function(i) {
  label_sid <- sample_order[[i]]
  value_sid <- value_order[[i]]
  data.table(
    Niche = niche_order,
    Sample = unname(sample_labels[[label_sid]]),
    rho = self_dt[[paste0(value_sid, "_rho")]],
    p = self_dt[[paste0(value_sid, "_p")]],
    FDR = self_dt[[paste0(value_sid, "_FDR")]]
  )
}))

plot_dt[, `:=`(
  Niche = factor(Niche, levels = niche_order),
  Sample = factor(Sample, levels = rev(unname(sample_labels[sample_order]))),
  sig = FDR < 0.05
)]

p <- ggplot(plot_dt, aes(Niche, Sample)) +
  geom_point(
    aes(size = -log10(pmax(FDR, 1e-300)), fill = rho, alpha = sig),
    shape = 21, color = "grey30", stroke = 0.25
  ) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
  scale_size(range = c(2, 6)) +
  guides(size = guide_legend(order = 1), fill = guide_colorbar(order = 2)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 10)
  ) +
  labs(fill = "Spearman rho", size = "-log10(FDR)", x = NULL, y = NULL)

ggsave(output_pdf, p, width = 5.6, height = 3.4,
       device = "pdf", useDingbats = FALSE)
