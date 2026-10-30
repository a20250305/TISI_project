suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source_csv <- "Figure_2H_source_data.csv"
output_pdf <- "Figure_2H_reproduced.pdf"

dt <- fread(source_csv)
required <- c("Sample", "Dataset", "CancerType", "Stage", "TISI_raw")
stopifnot(all(required %in% names(dt)))

dt[, TISI_raw := as.numeric(TISI_raw)]
stopifnot(
  nrow(dt) == 1380L,
  uniqueN(dt$Sample) == 1380L,
  uniqueN(dt$Dataset) == 15L,
  uniqueN(dt$CancerType) == 4L,
  all(is.finite(dt$TISI_raw)),
  all(dt$Stage %chin% c("Pre", "Post"))
)

dt[, Stage := factor(Stage, levels = c("Pre", "Post"))]

stat <- dt[, {
  pre <- TISI_raw[Stage == "Pre"]
  post <- TISI_raw[Stage == "Post"]
  p <- tryCatch(t.test(post, pre)$p.value, error = function(e) NA_real_)
  .(
    n_pre = length(pre),
    n_post = length(post),
    mean_pre = mean(pre),
    mean_post = mean(post),
    diff_mean = mean(post) - mean(pre),
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
labels <- unique(paste(dt$CancerType, dt$Stage, sep = "_"))
gaps <- paste0(unique(dt$CancerType), "_gap")
x_levels <- as.vector(rbind(
  labels[seq(1, length(labels), 2)],
  labels[seq(2, length(labels), 2)],
  gaps
))
x_levels <- x_levels[!is.na(x_levels)]

dt[, xlab := paste(CancerType, Stage, sep = "_")]
dt <- rbind(
  dt,
  data.table(
    CancerType = unique(dt$CancerType),
    Stage = "gap",
    xlab = gaps,
    TISI_raw = NA_real_
  ),
  fill = TRUE
)
dt[, xlab := factor(xlab, levels = x_levels)]

annotation <- merge(
  stat[, .(CancerType, sig)],
  dt[!is.na(TISI_raw), .(y = max(TISI_raw)), by = CancerType],
  by = "CancerType"
)
annotation[, `:=`(
  x = paste(CancerType, "Post", sep = "_"),
  y = y * 1.06
)]

plot <- ggplot(dt, aes(x = xlab, y = TISI_raw, fill = Stage)) +
  geom_violin(
    scale = "width", trim = FALSE, color = "black",
    linewidth = 0.25, na.rm = TRUE
  ) +
  geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA, na.rm = TRUE) +
  scale_fill_manual(
    values = c(Pre = "#3C5488", Post = "#DC0000"),
    na.translate = FALSE
  ) +
  geom_text(
    data = annotation,
    aes(x = x, y = y, label = sig),
    inherit.aes = FALSE,
    size = 3.6
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
      colour = "black", fill = NA, linewidth = 0.6
    ),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(
  filename = output_pdf,
  plot = plot,
  width = 7,
  height = 4,
  units = "in",
  device = "pdf"
)
