
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_file <- "Figure 2H_source_data.csv"
output_file <- "TISI_raw_violin_PrePost_byCancerType_reproduced_from_source.pdf"

dt <- fread(source_file)
required_cols <- c("Sample", "Dataset", "CancerType", "Stage", "TISI_raw")
stopifnot(all(required_cols %in% names(dt)))

dt[, TISI_raw := as.numeric(TISI_raw)]
stopifnot(
  nrow(dt) == 1380L,
  uniqueN(dt$Sample) == 1380L,
  uniqueN(dt$Dataset) == 15L,
  all(is.finite(dt$TISI_raw)),
  all(dt$Stage %chin% c("Pre", "Post")),
  uniqueN(dt$CancerType) == 4L
)

dt[, Stage := factor(Stage, levels = c("Pre", "Post"))]

## Same per-cancer-type Welch t-test and BH/FDR annotation as the original.
stat <- dt[, {
  x <- TISI_raw[Stage == "Pre"]
  y <- TISI_raw[Stage == "Post"]
  p <- tryCatch(t.test(y, x)$p.value, error = function(e) NA_real_)
  .(
    n_pre = length(x),
    n_post = length(y),
    mean_pre = mean(x),
    mean_post = mean(y),
    diff_mean = mean(y) - mean(x),
    p_mean = p
  )
}, by = CancerType]
stat[, FDR := p.adjust(p_mean, "BH")]
stat[, sig := fifelse(
  is.na(FDR), "NA",
  fifelse(FDR < 1e-3, "***",
          fifelse(FDR < 1e-2, "**",
                  fifelse(FDR < 0.05, "*", "ns")))
)]

setorder(dt, CancerType, Stage)
labs <- unique(paste(dt$CancerType, dt$Stage, sep = "_"))
gap <- paste0(unique(dt$CancerType), "_gap")
xlev <- as.vector(rbind(
  labs[seq(1, length(labs), 2)],
  labs[seq(2, length(labs), 2)],
  gap
))
xlev <- xlev[!is.na(xlev)]

dt[, xlab := paste(CancerType, Stage, sep = "_")]
dt <- rbind(
  dt,
  data.table(
    CancerType = unique(dt$CancerType),
    Stage = "gap",
    xlab = gap,
    TISI_raw = NA_real_
  ),
  fill = TRUE
)
dt[, xlab := factor(xlab, levels = xlev)]

ann <- merge(
  stat[, .(CancerType, sig)],
  dt[!is.na(TISI_raw), .(y = max(TISI_raw)), by = CancerType],
  by = "CancerType"
)
ann[, `:=`(
  x = paste(CancerType, "Post", sep = "_"),
  y = y * 1.06
)]

cols <- c(Pre = "#3C5488", Post = "#DC0000")

p <- ggplot(dt, aes(x = xlab, y = TISI_raw, fill = Stage)) +
  geom_violin(
    scale = "width", trim = FALSE, color = "black",
    linewidth = .25, na.rm = TRUE
  ) +
  geom_boxplot(width = .12, alpha = .9, outlier.shape = NA, na.rm = TRUE) +
  scale_fill_manual(values = cols, na.translate = FALSE) +
  geom_text(
    data = ann, aes(x = x, y = y, label = sig),
    inherit.aes = FALSE, size = 3.6
  ) +
  labs(
    x = NULL,
    y = "TISI (raw)",
    title = "TISI_raw across cancer types (Pre → Post)",
    subtitle = "Per-cancer-type mean comparison (Welch t-test, FDR-adjusted)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = .6
    ),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(output_file, p, width = 7, height = 4)
cat("Wrote:", output_file, "\n")
