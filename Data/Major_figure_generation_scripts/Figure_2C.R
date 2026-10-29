#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

source_csv <- "Figure_2C_source_data.csv"
output_pdf <- "Figure_2C_reproduced.pdf"

d <- fread(source_csv)
stopifnot(all(c("panel_label", "sample_id", "group", "TISI_minmax") %in% names(d)))

panel_order <- c(
  "GSE162175 - IMR90", "GSE284666 - BJ", "GSE301164 - BJ",
  "GSE305898 - iEC", "GSE305898 - HUVEC"
)
group_order <- c("Ctl", "TIS", "OIS", "RS")
group_labels <- c(Ctl = "Control", TIS = "TIS", OIS = "OIS", RS = "RS")
group_colors <- c(Ctl = "#0072B2", TIS = "#D62728", OIS = "#7B61A8", RS = "#B15928")

d[, panel_label := factor(panel_label, levels = panel_order)]
d[, panel_id := match(as.character(panel_label), panel_order)]
d[, group := factor(group, levels = group_order)]
d[, TISI_minmax := as.numeric(TISI_minmax)]

format_p <- function(p) {
  if (p < 0.001) sprintf("p=%.3e", p) else sprintf("p=%.4f", p)
}

pairwise <- rbindlist(lapply(panel_order, function(panel) {
  x <- d[panel_label == panel]
  groups <- group_order[group_order %in% unique(as.character(x$group))]
  rbindlist(lapply(seq_along(combn(groups, 2, simplify = FALSE)), function(i) {
    pair <- combn(groups, 2, simplify = FALSE)[[i]]
    a <- x[group == pair[1], TISI_minmax]
    b <- x[group == pair[2], TISI_minmax]
    p <- t.test(a, b, var.equal = TRUE)$p.value
    data.table(panel_label = panel, group1 = pair[1], group2 = pair[2],
               x1 = match(pair[1], groups), x2 = match(pair[2], groups),
               rank = i, p_label = format_p(p))
  }))
}))

pairwise[, panel_id := match(panel_label, panel_order)]
pairwise[, y := 1.05 + 0.10 * (rank - 1)]

p <- ggplot(d, aes(group, TISI_minmax, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.7) +
  geom_jitter(width = 0.10, height = 0, size = 1.8, alpha = 0.9) +
  geom_segment(data = pairwise,
               aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6, color = "black") +
  geom_segment(data = pairwise,
               aes(x = x1, xend = x1, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.6, color = "black") +
  geom_segment(data = pairwise,
               aes(x = x2, xend = x2, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.6, color = "black") +
  geom_text(data = pairwise,
            aes((x1 + x2) / 2, y + 0.015, label = p_label),
            inherit.aes = FALSE, size = 2.7, color = "black") +
  facet_wrap(
    ~panel_id,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(setNames(panel_order, as.character(seq_along(panel_order))))
  ) +
  scale_x_discrete(labels = group_labels, drop = TRUE) +
  scale_y_continuous(limits = c(0, 1.36), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "TISI Score") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(output_pdf, p, width = 12, height = 3.5, units = "in", device = "pdf")
