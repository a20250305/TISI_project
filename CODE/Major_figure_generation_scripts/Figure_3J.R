suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

source_csv <- "Figure_3J_source_data.csv"
output_pdf <- "KM_OS_TISI_scaled_reproduced_from_source_data_V2.pdf"

dt <- fread(source_csv)
stopifnot(all(c("patient_id", "TISI_scaled", "Censor", "OS") %in% names(dt)))

dt[, TISI_scaled := as.numeric(TISI_scaled)]
dt[, Censor := as.numeric(Censor)]
dt[, OS := as.numeric(OS)]
dt <- dt[is.finite(TISI_scaled) & is.finite(OS) & !is.na(Censor)]

cutoff <- median(dt$TISI_scaled)
dt[, grp := factor(
  ifelse(TISI_scaled >= cutoff, "TISI-high", "TISI-low"),
  levels = c("TISI-low", "TISI-high")
)]

fit <- survfit(Surv(OS, Censor) ~ grp, data = dt)
p_logrank <- 1 - pchisq(survdiff(Surv(OS, Censor) ~ grp, data = dt)$chisq, df = 1)

cox <- summary(coxph(Surv(OS, Censor) ~ grp, data = dt))
hr <- unname(cox$coef[1, "exp(coef)"])
ci_low <- unname(cox$conf.int[1, "lower .95"])
ci_high <- unname(cox$conf.int[1, "upper .95"])

n_low <- sum(dt$grp == "TISI-low")
n_high <- sum(dt$grp == "TISI-high")
events_low <- sum(dt$grp == "TISI-low" & dt$Censor == 1)
events_high <- sum(dt$grp == "TISI-high" & dt$Censor == 1)

subtitle <- sprintf(
  "Low n=%d (%d events) | High n=%d (%d events)\nHR=%.2f (%.2f-%.2f), log-rank P=%.3g",
  n_low, events_low, n_high, events_high,
  hr, ci_low, ci_high, p_logrank
)

km_plot <- ggsurvplot(
  fit,
  data = dt,
  conf.int = TRUE,
  conf.int.alpha = 0.25,
  palette = c("#2166AC", "#B2182B"),
  size = 1.2,
  xlab = "Overall survival (days)",
  ylab = "Survival probability",
  title = "GSE26253 OS ~ TISI_scaled",
  subtitle = subtitle,
  ggtheme = theme_classic(base_size = 16) +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.8
      ),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 14)
    )
)

ggsave(output_pdf, km_plot$plot, width = 5.5, height = 6, units = "in", device = "pdf")
