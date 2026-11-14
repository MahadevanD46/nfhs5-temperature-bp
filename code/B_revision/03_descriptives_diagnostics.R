# =============================================================================
# 03_descriptives_diagnostics.R
#   Table 1 with unweighted n and weighted %                    (comment 24)
#   State x year-month participant counts; phases                (comment 5)
#   Interview-to-biomarker gap distribution                      (comment 10)
#   Standardisation shift distribution and extrapolation         (comment 14)
#   Identifying variation after fixed effects                    (comment 15)
# =============================================================================
source("code/B_revision/00_config.R")
a <- load_a()
cc <- a[, complete.cases(t0, age_grp, hv104, hv270, hv025, shb14a, shb14b, shb14c, sh25, sh26)]
a <- a[cc]
W <- a$wt

## ---- Table 1 ---------------------------------------------------------------
cat_row <- function(label, x) {
  lev <- sort(unique(x[!is.na(x)]))
  rbindlist(lapply(lev, function(l) data.table(
    characteristic = label, level = as.character(l), n_unweighted = sum(x == l, na.rm = TRUE),
    pct_weighted = 100 * sum(W[x == l & !is.na(x)]) / sum(W[!is.na(x)]))))
}
num_row <- function(label, x) data.table(characteristic = label, level = "mean (SD)",
  n_unweighted = sum(!is.na(x)),
  pct_weighted = NA_real_, mean_weighted = wmean(x, W),
  sd_weighted = sqrt(Hmisc::wtd.var(x, W, na.rm = TRUE)))
t1 <- rbindlist(list(
  cat_row("Age group", a$age_grp), cat_row("Sex (1 men, 2 women)", a$hv104),
  cat_row("Residence (1 urban, 2 rural)", a$hv025), cat_row("Wealth quintile", a$hv270),
  cat_row("Current tobacco use", a$sh25), cat_row("Alcohol use", a$sh26),
  cat_row("Season of measurement", a$season), cat_row("Survey phase", a$shphase),
  cat_row("Measured elevated BP", a$htn_measured), cat_row("Elevated BP or medication", a$htn_report),
  num_row("Systolic BP, mmHg", a$sbp), num_row("Diastolic BP, mmHg", a$dbp),
  num_row("Random glucose, mg/dL", a$glu), num_row("Measurement-day temperature, degC", a$t0)),
  fill = TRUE)
fwrite(t1, file.path(OUT_DIR, "table1_v2.csv"))
reg("N_ANALYTIC", f0(nrow(a))); reg("N_CLUSTERS", f0(uniqueN(a$hv001)))
reg("N_DISTRICTS", f0(uniqueN(a$shdist))); reg("N_STATES", uniqueN(a$hv024))
reg("N_GLU_MISSING", f0(a[glu_med == 0 & is.na(glu), .N]))
reg("N_GLU_MED", f0(a[glu_med == 1, .N]))

## ---- (5) state x year-month ---------------------------------------------------
sm <- dcast(a[, .N, by = .(state = state_lab[as.character(hv024)], ym)], state ~ ym,
            value.var = "N", fill = 0)
sm[, total := rowSums(.SD), .SDcols = -1]
fwrite(sm, file.path(OUT_DIR, "state_by_month_counts.csv"))
both <- a[, .(before = any(date < as.Date("2020-04-01") & date >= as.Date("2020-01-01")),
              after = any(date >= as.Date("2020-11-01"))), by = hv024]
reg("N_STATES_BOTH_SIDES", both[before & after, .N])
reg("STATES_BOTH_SIDES", paste(state_lab[as.character(both[before & after, hv024])], collapse = ", "))
ph <- a[, .(n = .N, first = min(date), last = max(date)), by = shphase]
fwrite(ph, file.path(OUT_DIR, "phase_dates.csv"))
reg("DATE_FIRST", format(min(a$date), "%d %B %Y")); reg("DATE_LAST", format(max(a$date), "%d %B %Y"))

## ---- (10) gap distribution ------------------------------------------------------
gp <- a[, .(n = .N), by = gap_cat][order(gap_cat)][, pct := 100 * n / sum(n)]
gq <- a[, as.list(quantile(gap, c(0, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE))]
fwrite(gp, file.path(OUT_DIR, "gap_distribution.csv")); fwrite(gq, file.path(OUT_DIR, "gap_quantiles.csv"))
for (g in levels(a$gap_cat)) reg(paste0("GAP_", gsub("[^0-9A-Za-z]", "_", g)), f1(gp[gap_cat == g, pct]))
reg("GAP_NEG_N", f0(a[gap < 0, .N])); reg("GAP_MAX", a[, max(abs(gap), na.rm = TRUE)])
reg("GAP_P95", gq[["95%"]]); reg("SAME_DAY_PCT", f1(100 * mean(a$same_day)))

## ---- (14) shift distribution and extrapolation -------------------------------
s <- a[!is.na(t_annual)]
q <- Hmisc::wtd.quantile(s$shift, s$wt, c(.05, .25, .5, .75, .95))
reg("SHIFT_MEDIAN", f1(q[3])); reg("SHIFT_IQR", paste0(f1(q[2]), " to ", f1(q[4])))
reg("SHIFT_P5_P95", paste0(f1(q[1]), " to ", f1(q[5])))
reg("SHIFT_RANGE", paste0(f1(min(s$shift)), " to ", f1(max(s$shift))))
rng <- a[, .(tmin = min(t0), tmax = max(t0), p01 = quantile(t0, .01), p99 = quantile(t0, .99)), by = shdist]
s <- merge(s, rng, by = "shdist")
s[, `:=`(out_range = t_annual < tmin | t_annual > tmax, out_p = t_annual < p01 | t_annual > p99)]
ex <- s[, .(n = .N, w_mean_shift = wmean(shift, wt), pct_outside_range = 100 * mean(out_range),
            pct_outside_p1_p99 = 100 * mean(out_p)), by = .(state = state_lab[as.character(hv024)])]
fwrite(ex[order(-pct_outside_range)], file.path(OUT_DIR, "shift_extrapolation_by_state.csv"))
reg("EXTRAP_PCT", f1(100 * mean(s$out_range))); reg("EXTRAP_PCT_P", f1(100 * mean(s$out_p)))
reg("EXTRAP_STATES", paste(ex[pct_outside_range > 5][order(-pct_outside_range), state], collapse = ", "))
fwrite(data.table(shift = s$shift, wt = s$wt)[sample(.N, min(.N, 2e5))],
       file.path(OUT_DIR, "shift_sample_for_histogram.csv"))
rm(s); invisible(gc())

## ---- (15) identifying variation ----------------------------------------------
r1 <- feols(t0 ~ 1 | shdist + ym + dow, a, weights = ~wt, notes = FALSE)
r2 <- feols(t0 ~ 1 | shdist^ym + dow, a, weights = ~wt, notes = FALSE)
a[, t_res := resid(r1, na.rm = FALSE)]
dm <- a[, .(n = .N, dates = uniqueN(date), trange = max(t0) - min(t0)), by = .(shdist, ym)]
idv <- data.table(
  raw_sd = sqrt(Hmisc::wtd.var(a$t0, a$wt)),
  resid_sd_main_fe = sqrt(Hmisc::wtd.var(resid(r1), a$wt[obs(r1)])),
  resid_sd_dm_fe = sqrt(Hmisc::wtd.var(resid(r2), a$wt[obs(r2)])),
  district_months = nrow(dm),
  pct_dm_2plus_dates = 100 * mean(dm$dates >= 2),
  pct_obs_in_dm_2plus = 100 * dm[dates >= 2, sum(n)] / sum(dm$n),
  median_dates_per_dm = median(dm$dates),
  iqr_dates = paste(quantile(dm$dates, c(.25, .75)), collapse = "-"),
  median_within_dm_range = median(dm[dates >= 2, trange]),
  districts_with_within_dm_variation = dm[dates >= 2 & trange > 0, uniqueN(shdist)],
  districts_total = uniqueN(a$shdist))
fwrite(idv, file.path(OUT_DIR, "identifying_variation.csv")); print(idv)
dd <- dm[, .N, by = .(dates_cat = cut(dates, c(0, 1, 4, 9, 19, Inf), labels = c("1", "2-4", "5-9", "10-19", ">=20")))]
fwrite(dd[order(dates_cat)][, pct := 100 * N / sum(N)], file.path(OUT_DIR, "dates_per_district_month.csv"))
reg("ID_RAW_SD", f1(idv$raw_sd)); reg("ID_RES_SD", f2(idv$resid_sd_main_fe))
reg("ID_RES_SD_DM", f2(idv$resid_sd_dm_fe)); reg("ID_DM_N", f0(idv$district_months))
reg("ID_DM_2PLUS", f1(idv$pct_dm_2plus_dates)); reg("ID_OBS_2PLUS", f1(idv$pct_obs_in_dm_2plus))
reg("ID_MED_DATES", idv$median_dates_per_dm); reg("ID_IQR_DATES", idv$iqr_dates)
reg("ID_MED_RANGE", f1(idv$median_within_dm_range))
reg("ID_DIST_VAR", f0(idv$districts_with_within_dm_variation))
reg("ID_DIST_TOT", f0(idv$districts_total))
saveRDS(a[, .(hv001, hv002, hvidx, t_res)], file.path(OUT_DIR, "t_residuals.rds"))
log_msg("03 done")
