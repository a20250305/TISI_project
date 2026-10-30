source_csv <- "Figure_1A_source_data.csv"
output_pdf <- "Figure_1A_from_source_data.pdf"

dt <- read.csv(source_csv, check.names = FALSE)
stopifnot(all(c("cancer_type", "total_n", "control_n", "sene_n") %in% names(dt)))
stopifnot(all(dt$total_n == dt$control_n + dt$sene_n))
stopifnot(sum(dt$total_n) == 258, sum(dt$control_n) == 106, sum(dt$sene_n) == 152)

# Figure colors are defined in the code, not stored in the source data.
cancer_colors <- c(
  "Colorectal cancer" = "#526D8E",
  "DLBCL" = "#C76F52",
  "Glioma" = "#5A8F72",
  "HR+ Breast cancer" = "#A6546F",
  "Liposarcoma" = "#C5A142",
  "Liver cancer" = "#4F918B",
  "Lung cancer" = "#6E6AA8",
  "Melanoma" = "#D08A8A",
  "Osteosarcoma" = "#8F6A57",
  "Pancreatic cancer" = "#D1773A",
  "Angiosarcoma" = "#6A8747",
  "TN Breast cancer" = "#956084"
)
cancer_fill <- unname(cancer_colors[dt$cancer_type])
stopifnot(all(!is.na(cancer_fill)))

control_color <- "#3B6FB6"
sene_color <- "#D55E5E"

sector_polygon <- function(cx, cy, r0, r1, a0, a1, fill) {
  a <- seq(a0, a1, length.out = 80) * pi / 180
  outer <- cbind(cx + r1 * cos(a), cy + r1 * sin(a))
  inner <- cbind(cx + r0 * cos(rev(a)), cy + r0 * sin(rev(a)))
  polygon(rbind(outer, inner), col = fill, border = "white", lwd = 0.9)
}

draw_legend <- function(x, y_top, labels, colors) {
  for (i in seq_along(labels)) {
    y <- y_top - (i - 1) * 0.285
    rect(x, y - 0.09, x + 0.20, y + 0.09, col = colors[i], border = NA)
    text(x + 0.30, y, labels[i], adj = c(0, 0.5), cex = 0.70, col = "#222222")
  }
}

pdf(output_pdf, width = 8.4, height = 5.0, family = "Helvetica", useDingbats = FALSE)
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bty = "n", xpd = NA, family = "Helvetica")
plot.new()
plot.window(xlim = c(0, 8.4), ylim = c(0, 5.0), asp = 1)
text(0.18, 4.55, "A", adj = c(0, 0.5), font = 2, cex = 1.35, col = "#222222")

outer_order <- rev(seq_len(nrow(dt)))
outer_n <- dt$total_n[outer_order]
angles <- 90 - c(0, cumsum(360 * outer_n / sum(outer_n)))
cx <- 1.82
cy <- 2.88
r_outer <- 1.27
r_inner <- 0.78
r_hole <- 0.48

for (i in seq_along(outer_order)) {
  sector_polygon(cx, cy, r_inner, r_outer, angles[i], angles[i + 1], cancer_fill[outer_order[i]])
}

for (i in seq_along(outer_order)) {
  j <- outer_order[i]
  split <- angles[i] - (angles[i] - angles[i + 1]) * dt$control_n[j] / dt$total_n[j]
  sector_polygon(cx, cy, r_hole, r_inner, angles[i], split, control_color)
  sector_polygon(cx, cy, r_hole, r_inner, split, angles[i + 1], sene_color)
}

theta <- seq(0, 2 * pi, length.out = 200)
polygon(cx + r_hole * cos(theta), cy + r_hole * sin(theta), col = "white", border = NA)
text(cx, cy + 0.12, format(sum(dt$total_n), big.mark = ","), cex = 1.00, col = "#222222")
text(cx, cy - 0.17, "samples", cex = 1.00, col = "#222222")

text(4.03, 4.25, "Cancer type", adj = c(0, 0.5), cex = 0.82, col = "#222222")
draw_legend(4.03, 3.98, dt$cancer_type, cancer_fill)

legend_y <- 0.58
rect(1.03, legend_y - 0.09, 1.23, legend_y + 0.09, col = control_color, border = NA)
text(1.34, legend_y, "Control", adj = c(0, 0.5), cex = 0.72, col = "#222222")
rect(2.12, legend_y - 0.09, 2.32, legend_y + 0.09, col = sene_color, border = NA)
text(2.43, legend_y, "Sene", adj = c(0, 0.5), cex = 0.72, col = "#222222")

dev.off()
