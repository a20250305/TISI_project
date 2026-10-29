source_csv <- "Figure_1C_source_data.csv"
output_pdf <- "Figure_1C_from_source_data.pdf"

dt <- read.csv(source_csv, check.names = FALSE)
required <- c("sample_id", "cancer_type", "agent_category", "cell_line")
stopifnot(all(required %in% names(dt)))
stopifnot(!anyDuplicated(dt$sample_id))

cancer_order <- c(
  "Colorectal cancer", "DLBCL", "Glioma", "HR+ Breast cancer",
  "Liposarcoma", "Liver cancer", "Lung cancer", "Melanoma",
  "Osteosarcoma", "Pancreatic cancer", "Angiosarcoma", "TN Breast cancer"
)
agent_order <- c(
  "Radiation", "Mitotic kinase inhibitor", "Chemotherapy", "Cell cycle inhibitor"
)
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
stopifnot(all(dt$cancer_type %in% cancer_order))
stopifnot(all(dt$agent_category %in% agent_order))
stopifnot(all(cancer_order %in% names(cancer_colors)))

counts <- aggregate(
  sample_id ~ cancer_type + agent_category,
  data = dt,
  FUN = length
)
names(counts)[3] <- "n"
counts$cancer_type <- factor(counts$cancer_type, levels = cancer_order)
counts$agent_category <- factor(counts$agent_category, levels = agent_order)
counts <- counts[order(counts$agent_category, counts$cancer_type), ]

counts$fill <- unname(cancer_colors[as.character(counts$cancer_type)])
stopifnot(!anyNA(counts$fill), nrow(counts) == 22, sum(counts$n) == 150)

size_n <- c(3, 4, 5, 6, 7, 8, 9, 10, 12, 18, 24)
size_cex <- c(0.30, 0.62, 0.75, 0.84, 0.93, 1.01,
              1.08, 1.17, 1.24, 1.50, 1.53)
cex_for_count <- function(n) {
  approx(size_n, size_cex, xout = n, rule = 2)$y
}

pdf(output_pdf, width = 270 / 72, height = 110 / 72,
    family = "Helvetica", useDingbats = FALSE)
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bty = "n", xpd = NA,
    family = "Helvetica")
plot.new()
plot.window(xlim = c(0, 270), ylim = c(0, 110))

plot_x0 <- 88
plot_x1 <- 236
plot_y0 <- 50
plot_y1 <- 97
x_pos <- seq(plot_x0 + 7, plot_x1 - 7, length.out = length(cancer_order))
y_pos <- setNames(seq(90, 60, length.out = length(agent_order)), agent_order)

text(5, 103, "C", adj = c(0, 0.5), cex = 0.80, font = 2, col = "#20252B")
abline(h = unname(y_pos), v = x_pos,
       col = adjustcolor("#D8DDE2", alpha.f = 0.68), lwd = 0.55)
rect(plot_x0, plot_y0, plot_x1, plot_y1, border = "#5C636B", lwd = 0.65)

counts$x <- x_pos[match(as.character(counts$cancer_type), cancer_order)]
counts$y <- unname(y_pos[as.character(counts$agent_category)])
points(counts$x, counts$y, pch = 16,
       cex = cex_for_count(counts$n), col = counts$fill)

text(plot_x0 - 2, unname(y_pos), names(y_pos),
     adj = c(1, 0.5), cex = 0.55, col = "#30353A")
text(x_pos, 47, cancer_order, srt = 45, adj = c(1, 1),
     cex = 0.53, col = "#30353A")

legend_x <- 251
legend_y <- c(88, 79, 70, 61, 52)
legend_n <- c(3, 6, 9, 12, 18)
text(legend_x, 98, "Count", adj = c(0.5, 0.5),
     cex = 0.55, col = "#353A40")
points(rep(legend_x, length(legend_n)), legend_y, pch = 16,
       cex = cex_for_count(legend_n), col = "#333333")
text(legend_x + 7, legend_y, legend_n,
     adj = c(0, 0.5), cex = 0.50, col = "#30353A")

dev.off()
