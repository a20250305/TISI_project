

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  dirname(normalizePath(script_path))
} else {
  getwd()
}

source_file <- file.path(
  script_dir, "..", "Figure 2H",
  "Figure 2H_source_data_TISI_raw_violin_PrePost_byCancerType.csv"
)
output_file <- file.path(
  script_dir,
  "Figure_2F_CancerType_donut_outer_ring_original_Donut_style.pdf"
)

dt <- fread(source_file)
required_cols <- c("Sample", "Dataset", "CancerType", "Stage", "TISI_raw")
stopifnot(all(required_cols %in% names(dt)))

dt[, TISI_raw := as.numeric(TISI_raw)]
dt <- dt[
  is.finite(TISI_raw) &
    Stage %chin% c("Pre", "Post")
]

## Audit checks for the released source data.
stopifnot(
  nrow(dt) == 1380L,
  uniqueN(dt$Sample) == 1380L,
  uniqueN(dt$Dataset) == 15L,
  uniqueN(dt$CancerType) == 4L
)

## These lines intentionally match Code_6.1_clinical_validation_022226.R.
ct_tab <- dt[, .(n = .N), by = CancerType][order(CancerType)]
cancers <- ct_tab$CancerType
basecol <- setNames(
  viridis(length(cancers), option = "D", begin = .05, end = .85),
  cancers
)

ct_tab[, ymax := cumsum(n)]
ct_tab[, ymin := shift(ymax, fill = 0)]

p <- ggplot() +
  ## Only the original CancerType outer ring is retained.
  geom_rect(
    data = ct_tab,
    aes(
      ymin = ymin, ymax = ymax,
      xmin = 2.2, xmax = 3,
      fill = CancerType
    ),
    color = "white",
    linewidth = .8
  ) +
  coord_polar(theta = "y") +
  xlim(.2, 3) +
  scale_fill_manual(values = basecol, breaks = cancers) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  filename = output_file,
  plot = p,
  width = 8,
  height = 7,
  units = "in"
)

cat("DONE ->", output_file, "\n")
print(ct_tab[, .(CancerType, n)])
