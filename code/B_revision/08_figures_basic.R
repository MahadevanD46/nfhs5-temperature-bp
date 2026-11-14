# =============================================================================
# 08_figures_basic.R — basic plots of the main results (default ggplot styling)
# Content corresponds to Figures 2-3 and Supplementary Figures S1-S5.
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages(library(ggplot2))
sv <- function(g, f) ggsave(file.path(FIG_DIR, f), g)

# spline curves (all outcomes) with category estimates for SBP
sc <- fread(file.path(OUT_DIR, "spline_curves.csv"))
sv(ggplot(sc, aes(temp, est)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) + geom_line() +
     facet_wrap(~outcome, scales = "free_y"), "splines.png")
ct <- fread(file.path(OUT_DIR, "temperature_categories.csv"))
sv(ggplot(ct[outcome == "sbp"], aes(mean_temp, est)) + geom_pointrange(aes(ymin = lo, ymax = hi)),
   "sbp_categories.png")

# effect modification: season slopes
em <- fread(file.path(OUT_DIR, "effect_modification_slopes.csv"))
sv(ggplot(em[modifier == "Season" & outcome == "sbp"], aes(level, est)) +
     geom_pointrange(aes(ymin = lo, ymax = hi)), "season_slopes.png")

# fieldwork calendar
a <- load_a(cols = c("date", "shphase"))
sv(ggplot(a[, .N, by = .(month = as.Date(format(date, "%Y-%m-01")), shphase)],
          aes(month, N, fill = factor(shphase))) + geom_col(), "calendar.png")
rm(a)

# standardisation shift
sh <- fread(file.path(OUT_DIR, "shift_sample_for_histogram.csv"))
sv(ggplot(sh, aes(shift, weight = wt)) + geom_histogram(binwidth = 1), "shift.png")

# state-specific slopes
ss <- fread(file.path(OUT_DIR, "state_slopes.csv"))
sv(ggplot(ss, aes(est, state)) + geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se)) +
     facet_wrap(~outcome), "state_slopes.png")
log_msg("08 done")
