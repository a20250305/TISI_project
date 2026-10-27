source_csv <- "Figure_1B_right_source_data.csv"
output_pdf <- "Figure_1B_right_from_source_data.pdf"

dt <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_id", "cell_line", "cancer_type", "dataset")
stopifnot(all(required %in% names(dt)))
stopifnot(!anyDuplicated(dt$sample_id))

cancer_order <- c(
  "Colorectal cancer", "DLBCL", "Glioma", "HR+ Breast cancer",
  "Liposarcoma", "Liver cancer", "Lung cancer", "Melanoma",
  "Osteosarcoma", "Pancreatic cancer", "Angiosarcoma", "TN Breast cancer"
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

stopifnot(nrow(dt) == 258)
stopifnot(length(unique(dt$cell_line)) == 40)
stopifnot(length(unique(dt$cancer_type)) == 12)
stopifnot(length(unique(dt$dataset)) == 14)
stopifnot(all(dt$cancer_type %in% cancer_order))

assoc <- unique(dt[c("cell_line", "cancer_type", "dataset")])
dataset_order <- unique(assoc$dataset)
assoc$cancer_idx <- match(assoc$cancer_type, cancer_order)
assoc$dataset_idx <- match(assoc$dataset, dataset_order)
assoc <- assoc[order(assoc$cancer_idx, assoc$cell_line, assoc$dataset_idx), ]

cancer_counts <- aggregate(
  sample_id ~ cancer_type,
  data = dt,
  FUN = length
)
names(cancer_counts)[2] <- "n_samples"
cancer_counts$cancer_type <- factor(
  cancer_counts$cancer_type,
  levels = cancer_order
)
cancer_counts <- cancer_counts[order(cancer_counts$cancer_type), ]

dataset_counts <- aggregate(
  sample_id ~ dataset,
  data = dt,
  FUN = length
)
names(dataset_counts)[2] <- "n_samples"
dataset_counts$dataset <- factor(
  dataset_counts$dataset,
  levels = dataset_order
)
dataset_counts <- dataset_counts[order(dataset_counts$dataset), ]

stopifnot(sum(cancer_counts$n_samples) == 258)
stopifnot(!anyNA(cancer_counts$n_samples), !anyNA(dataset_counts$n_samples))

x_mid0 <- 0.15
x_mid1 <- 1.30
x_right0 <- 3.30
x_right1 <- 4.22
y_bottom <- 0.95
y_top <- 5.65

mid_gap <- 0.015
mid_scale <- (y_top - y_bottom - mid_gap * (nrow(cancer_counts) - 1)) /
  sum(cancer_counts$n_samples)
cancer_counts$y0 <- NA_real_
cancer_counts$y1 <- NA_real_
cancer_counts$y <- NA_real_
cursor <- y_top
for (i in seq_len(nrow(cancer_counts))) {
  h <- cancer_counts$n_samples[i] * mid_scale
  cancer_counts$y1[i] <- cursor
  cancer_counts$y0[i] <- cursor - h
  cancer_counts$y[i] <- (cancer_counts$y0[i] + cancer_counts$y1[i]) / 2
  cursor <- cursor - h - mid_gap
}

dataset_top_gap <- 0.27
dataset_top <- y_top - dataset_top_gap
dataset_scale <- (dataset_top - y_bottom) / sum(dataset_counts$n_samples)
right_heights <- dataset_counts$n_samples * dataset_scale
dataset_counts$y0 <- NA_real_
dataset_counts$y1 <- NA_real_
dataset_counts$y <- NA_real_
cursor <- dataset_top
for (i in seq_len(nrow(dataset_counts))) {
  h <- right_heights[i]
  dataset_counts$y1[i] <- cursor
  dataset_counts$y0[i] <- cursor - h
  dataset_counts$y[i] <- (dataset_counts$y0[i] + dataset_counts$y1[i]) / 2
  cursor <- cursor - h
}

assoc$mid_y <- NA_real_
assoc$right_y <- NA_real_
assoc$mid_h <- NA_real_
assoc$right_h <- NA_real_

for (ct in cancer_order) {
  rows <- which(assoc$cancer_type == ct)
  block <- cancer_counts[cancer_counts$cancer_type == ct, ]
  spacing <- (block$y1 - block$y0) / (length(rows) + 1)
  assoc$mid_y[rows] <- block$y0 + rev(seq_along(rows)) * spacing
  assoc$mid_h[rows] <- min(0.045, (block$y1 - block$y0) * 0.72 / length(rows))
}

for (ds in dataset_order) {
  rows <- which(assoc$dataset == ds)
  block <- dataset_counts[dataset_counts$dataset == ds, ]
  spacing <- (block$y1 - block$y0) / (length(rows) + 1)
  assoc$right_y[rows] <- block$y0 + rev(seq_along(rows)) * spacing
  assoc$right_h[rows] <- min(0.040, (block$y1 - block$y0) * 0.68 / length(rows))
}

bezier_curve <- function(x0, y0, x1, y1, n = 36) {
  t <- seq(0, 1, length.out = n)
  c1x <- x0 + (x1 - x0) * 0.35
  c2x <- x0 + (x1 - x0) * 0.65
  x <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * c1x +
    3 * (1 - t) * t^2 * c2x + t^3 * x1
  y <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * y0 +
    3 * (1 - t) * t^2 * y1 + t^3 * y1
  cbind(x, y)
}

ribbon <- function(x0, y0, x1, y1, h0, h1) {
  top <- bezier_curve(x0, y0 + h0 / 2, x1, y1 + h1 / 2)
  bottom <- bezier_curve(x0, y0 - h0 / 2, x1, y1 - h1 / 2)
  polygon(
    rbind(top, bottom[nrow(bottom):1, ]),
    col = "#F1F1F2",
    border = "#231F20",
    lwd = 0.45
  )
}

text_colour_for_fill <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  luminance <- 0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]
  if (luminance > 0.63) "#222222" else "white"
}

pdf(output_pdf, width = 3.25, height = 3.35,
    family = "Helvetica", useDingbats = FALSE)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0),
    xaxs = "i", yaxs = "i", bty = "n", xpd = NA,
    family = "Helvetica")
plot.new()
plot.window(xlim = c(-0.35, 4.45), ylim = c(0, 6.6))

for (i in seq_len(nrow(assoc))) {
  ribbon(x_mid1, assoc$mid_y[i], x_right0, assoc$right_y[i],
         assoc$mid_h[i], assoc$right_h[i])
}

for (i in seq_len(nrow(cancer_counts))) {
  ct <- as.character(cancer_counts$cancer_type[i])
  block_height <- cancer_counts$y1[i] - cancer_counts$y0[i]
  label_cex <- min(0.52, max(0.24, block_height * 1.35))
  rect(
    x_mid0, cancer_counts$y0[i], x_mid1, cancer_counts$y1[i],
    col = cancer_colors[ct], border = "#FFFFFF", lwd = 0.75
  )
  text(
    (x_mid0 + x_mid1) / 2, cancer_counts$y[i], ct,
    cex = label_cex,
    col = text_colour_for_fill(cancer_colors[ct])
  )
}

for (i in seq_len(nrow(dataset_counts))) {
  dataset_height <- dataset_counts$y1[i] - dataset_counts$y0[i]
  rect(
    x_right0, dataset_counts$y0[i], x_right1, dataset_counts$y1[i],
    col = "#FAFAFA", border = "#666A70", lwd = 0.55
  )
  text(
    (x_right0 + x_right1) / 2, dataset_counts$y[i],
    as.character(dataset_counts$dataset[i]),
    cex = min(0.44, max(0.24, dataset_height * 2.4)),
    col = "#333333"
  )
}

text(-0.30, 6.25, "B", adj = c(0, 0.5), font = 2, cex = 1.20,
     col = "#222222")
text((x_mid0 + x_mid1) / 2, 6.12, "Cancer type (12)",
     cex = 0.80, col = "#222222")
text((x_right0 + x_right1) / 2, 6.12, "Dataset (14)",
     cex = 0.80, col = "#222222")

dev.off()
