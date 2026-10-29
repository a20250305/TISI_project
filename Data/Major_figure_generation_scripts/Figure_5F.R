suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source_csv <- "Figure_5F_source_data.csv"
output_pdf <- "Figure_5F.pdf"

dt <- fread(source_csv)
stopifnot(all(c(
  "panel", "variable", "displayed_row_label", "rho",
  "FDR", "neglog10_FDR", "sig"
) %in% names(dt)))

sample_order <- c("MRD-1", "MRD-2", "MRD-3", "MRD-4")
niche_order <- paste0("Niche_", 1:13)
celltype_order <- c(
  "Tumor_cells", "Macrophages", "T_cells", "B_cells",
  "Dendritic_cells", "Endothelial_cells", "Fibroblasts"
)

dt[, displayed_row_label := factor(displayed_row_label, levels = rev(sample_order))]
dt[, sig := as.logical(sig)]

f5f_theme <- theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

rho_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
  limits = c(-0.3, 0.1), oob = scales::squish,
  na.value = "white", name = "Spearman rho",
  breaks = c(-0.3, -0.2, -0.1, 0, 0.1)
)

size_scale <- scale_size(
  range = c(1.5, 6), breaks = c(1, 2, 3), name = "-log10(FDR)"
)

make_panel <- function(panel_name, x_order) {
  z <- dt[panel == panel_name]
  z[, variable := factor(variable, levels = x_order)]
  ggplot(z, aes(x = variable, y = displayed_row_label)) +
    geom_point(
      aes(size = neglog10_FDR, fill = rho, alpha = sig),
      shape = 21, color = "grey30", stroke = 0.25, na.rm = TRUE
    ) +
    rho_scale +
    size_scale +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    f5f_theme
}

p_niche <- make_panel("Neighbor niche", niche_order)
p_celltype <- make_panel("Neighbor major cell type", celltype_order)

figure <- (p_niche | p_celltype) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "TISI associated tumor neighbor niche and cell types") &
  theme(
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5,
                              margin = margin(b = 8))
  )

ggsave(
  output_pdf,
  figure,
  width = 8.5,
  height = 3.8,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)
