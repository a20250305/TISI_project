suppressPackageStartupMessages(library(ggplot2))

source_csv <- "Figure_3H_source_data.csv"
output_pdf <- "Figure_3H_reproduced.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE)
stopifnot(all(c(
  "cancer_type", "selected_TISI_variable", "OR",
  "CI_low", "CI_high", "p"
) %in% names(dt)))
stopifnot(nrow(dt) == 2, all(is.finite(dt$OR)))
stopifnot(
  dt$selected_TISI_variable[match("COAD", dt$cancer_type)] == "TISI_scaled",
  dt$selected_TISI_variable[match("OV", dt$cancer_type)] == "TISI_raw"
)

dt$cancer_type <- factor(dt$cancer_type, levels = c("OV", "COAD"))
dt$label <- sprintf(
  "OR=%.2f\nP=%s",
  dt$OR,
  formatC(dt$p, format = "fg", digits = 3)
)

p <- ggplot(dt, aes(x = OR, y = cancer_type)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.8
  ) +
  geom_point(
    shape = 21,
    size = 4.2,
    fill = "#B2182B",
    colour = "black",
    stroke = 0.3
  ) +
  geom_text(
    aes(x = 1.08, y = cancer_type, label = label),
    hjust = 0,
    lineheight = 1.05,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2),
    labels = c("0.25", "0.50", "1.00", "2.00")
  ) +
  coord_cartesian(xlim = c(0.18, 2.2), clip = "off") +
  labs(
    x = "Odds Ratio",
    y = NULL,
    title = "High TISI predicts poor response (OR < 1)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.6
    ),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.margin = margin(6, 70, 6, 6)
  )

ggsave(
  output_pdf,
  p,
  width = 5.2,
  height = 2.5,
  units = "in",
  device = "pdf"
)
