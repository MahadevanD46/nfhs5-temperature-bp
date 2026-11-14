# =============================================================================
# 01_flow_validation.R — participant flow, excluded-vs-included comparison,
#                        validation of the BP derivation, biomarker-date handling
# Reviewer comments: 7 (flow, differential exclusion), 8 (validation), 10 (dates)
# Inputs : IAPR7EFL.DTA, IAGE7AFL.shp, output_full/analysis_dataset.rds
# Outputs: output_revision/flow_counts.csv            -> Supplementary Figure S2
#          output_revision/flow_by_state.csv          -> Supplementary Table S3b
#          output_revision/excluded_vs_included.csv   -> Supplementary Table S3
#          output_revision/validation_bp.csv          -> Supplementary Table S4
#          output_revision/date_handling.csv          -> Supplementary Table S5 (part)
# Memory-lean version: one copy of the PR file, steps marked with flags.
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages({ library(haven); library(lubridate); library(sf) })

## ---- 1. analytic sample keys first (then free the large analysis file) --------
a <- readRDS(file.path(IN_DIR, "analysis_dataset.rds")); setDT(a)
n_temp <- nrow(a)                                         # temperature assigned
cc <- complete.cases(a[, .(t0, shdist, ym, dow, wt, age_grp, hv104, hv270, hv025,
                           shb14a, shb14b, shb14c, sh25, sh26)])
a_keys <- a[cc, paste(hv001, hv002, hvidx)]
n_analytic <- length(a_keys)
rm(a); invisible(gc())

## ---- 2. PR file (only the columns needed) -------------------------------------
keep <- c("hv001","hv002","hvidx","hv005","hv103","hv006","hv007","hv016","hv024","hv025",
          "hv104","hv105","hv270","hv807d","hv807m","hv807y",
          "shb18s","shb18d","shb21","shb25s","shb25d","shb29s","shb29d","ha54")
pr <- read_dta(PR_FILE, col_select = tidyselect::all_of(keep)) |> zap_labels() |> as.data.table()
ok <- function(x, lo, hi) fifelse(!is.na(x) & x >= lo & x <= hi, as.numeric(x), NA_real_)

# identical derivations to the original pipeline
pr[, `:=`(s2 = ok(shb25s, 70, 270), s3 = ok(shb29s, 70, 270),
          d2 = ok(shb25d, 40, 150), d3 = ok(shb29d, 40, 150),
          bp_med = as.integer(fcoalesce(as.numeric(shb21), 0) == 1),
          preg = as.integer(hv104 == 2 & fcoalesce(as.numeric(ha54), 0) == 1),
          wt = hv005 / 1e6)]
pr[, sbp := rowMeans(cbind(s2, s3), na.rm = TRUE)][is.nan(sbp), sbp := NA]
pr[, dbp := rowMeans(cbind(d2, d3), na.rm = TRUE)][is.nan(dbp), dbp := NA]
pr[, bp_any_raw := !is.na(shb25s) | !is.na(shb29s) | !is.na(shb25d) | !is.na(shb29d)]
pr[, bio_date := make_date(ok(hv807y, 2019, 2021), ok(hv807m, 1, 12), ok(hv807d, 1, 31))]
pr[, int_date := make_date(hv007, hv006, hv016)]
pr[, date := fcoalesce(bio_date, int_date)]
pr[, htn_report := as.integer(sbp >= 140 | dbp >= 90 | bp_med == 1)]
pr[, htn_measured := as.integer(sbp >= 140 | dbp >= 90)]
pr[, valid_bp := !is.na(sbp) & !is.na(dbp)]

## ============================================================ (8) validation
# NFHS-5 report definition: age >=15, SBP >=140 or DBP >=90 or taking antihypertensive
# medication; BP = mean of readings 2 and 3; weights hv005.
PUBLISHED <- c(men = 24.0, women = 21.3)   # NFHS-5 India report (IIPS & ICF 2021)
pr[, s1 := ok(shb18s, 70, 270)][, d1 := ok(shb18d, 40, 150)]
v15 <- pr[hv105 >= 15 & preg == 0 & valid_bp == TRUE,
          .(sex = fifelse(hv104 == 1, "men", "women"), wt, hv103, bp_med, htn_report, htn_measured,
            r123 = as.integer(rowMeans(cbind(s1, s2, s3), na.rm = TRUE) >= 140 |
                              rowMeans(cbind(d1, d2, d3), na.rm = TRUE) >= 90 | bp_med == 1))]
val <- v15[, .(
  denominator_unweighted = .N, weighted_n = sum(wt),
  report_def_r23 = 100 * weighted.mean(htn_report, wt),
  report_def_r123 = 100 * weighted.mean(r123, wt),
  measured_only_r23 = 100 * weighted.mean(htn_measured, wt),
  on_medication = 100 * weighted.mean(bp_med, wt),
  report_def_defacto = 100 * weighted.mean(htn_report[hv103 == 1], wt[hv103 == 1])), by = sex]
rm(v15); pr[, c("s1", "d1", "shb18s", "shb18d") := NULL]; invisible(gc())
val[, published := PUBLISHED[sex]]
val18 <- pr[hv105 >= 18 & preg == 0 & valid_bp == TRUE,
            .(adults18_report_def = 100 * weighted.mean(htn_report, wt),
              adults18_measured_only = 100 * weighted.mean(htn_measured, wt), n18 = .N),
            by = .(sex = fifelse(hv104 == 1, "men", "women"))]
val <- merge(val, val18, by = "sex")
fwrite(val, file.path(OUT_DIR, "validation_bp.csv")); print(val)
for (s in c("men", "women")) {
  reg(paste0("VAL_", toupper(s)), f1(val[sex == s, report_def_r23]))
  reg(paste0("VAL_N_", toupper(s)), f0(val[sex == s, denominator_unweighted]))
  reg(paste0("VAL_R123_", toupper(s)), f1(val[sex == s, report_def_r123]))
  reg(paste0("VAL_DEFACTO_", toupper(s)), f1(val[sex == s, report_def_defacto]))
}

## ============================================================ (7) flow (flags, no copies)
ge <- st_read(GE_FILE, quiet = TRUE) |> st_drop_geometry() |> as.data.table()
gps_ok <- ge[!(LATNUM == 0 & LONGNUM == 0), DHSCLUST]; rm(ge)
pr[, `:=`(f_age = hv105 >= 18)]
pr[, f_preg := f_age & preg == 0]
pr[, f_bp := f_preg & valid_bp]
pr[, f_date := f_bp & !format(date, "%Y_%m") %in% c("2020_04", "2020_10")]
pr[, f_gps := f_date & hv001 %in% gps_ok]
pr[, k := paste(hv001, hv002, hvidx)]
pr[, f_analytic := f_preg & k %in% a_keys]
rm(a_keys); invisible(gc())

n_nobp  <- pr[f_preg & !valid_bp & !bp_any_raw, .N]
n_badbp <- pr[f_preg & !valid_bp & bp_any_raw, .N]
flow <- data.table(
  step = c("Household members in PR file", "Aged >=18 years", "Not pregnant",
           "Valid systolic and diastolic BP", "Measurement date outside COVID-19 suspension",
           "Cluster has GPS coordinates",
           "Temperature assigned (ERA5-Land land pixel or coastal infill)",
           "Complete covariates (analytic sample)"),
  n = c(nrow(pr), pr[, sum(f_age)], pr[, sum(f_preg)], pr[, sum(f_bp)], pr[, sum(f_date)],
        pr[, sum(f_gps)], n_temp, n_analytic))
flow[, excluded_at_step := shift(n) - n]
flow[, detail := ""]
flow[4, detail := sprintf("BP not measured: %s; readings outside plausible range: %s", f0(n_nobp), f0(n_badbp))]
fwrite(flow, file.path(OUT_DIR, "flow_counts.csv")); print(flow)
for (i in seq_len(nrow(flow))) reg(sprintf("FLOW_N%d", i), f0(flow$n[i]))
for (i in 2:nrow(flow)) reg(sprintf("FLOW_EX%d", i), f0(flow$excluded_at_step[i]))
reg("FLOW_NOBP", f0(n_nobp)); reg("FLOW_BADBP", f0(n_badbp))

# by state
st <- pr[f_preg == TRUE, .(eligible = .N, valid_bp = sum(f_bp), with_gps = sum(f_gps),
                           analytic = sum(f_analytic)), by = hv024]
st[, `:=`(state = state_lab[as.character(hv024)],
          pct_valid_bp = 100 * valid_bp / eligible,
          pct_analytic_of_valid_bp = 100 * analytic / pmax(valid_bp, 1))]
fwrite(st[order(state)], file.path(OUT_DIR, "flow_by_state.csv"))

## ============================================================ excluded vs included
pr[f_preg == TRUE, group := fcase(!valid_bp, "No valid BP", !f_analytic, "Valid BP, excluded later",
                                  default = "Analytic sample")]
pr[, age_g := cut(hv105, c(17, 29, 39, 49, 59, 69, Inf),
                  labels = c("18-29", "30-39", "40-49", "50-59", "60-69", ">=70"))]
pr[, region := region_of[as.character(hv024)]]
desc <- function(var, lev) {
  out <- pr[f_preg == TRUE, .(pct_unw = 100 * mean(get(var) == lev, na.rm = TRUE),
                              pct_w = 100 * weighted.mean(get(var) == lev, wt, na.rm = TRUE)), by = group]
  out[, characteristic := paste0(var, " = ", lev)]
}
rows <- rbindlist(c(
  list(pr[f_preg == TRUE, .(characteristic = "n", pct_unw = as.numeric(.N), pct_w = NA_real_), by = group]),
  list(desc("hv104", 2)), list(desc("hv025", 1)),
  lapply(levels(pr$age_g), function(l) desc("age_g", l)),
  lapply(1:5, function(q) desc("hv270", q)),
  lapply(sort(unique(na.omit(pr$region))), function(r) desc("region", r))), fill = TRUE)
tab <- dcast(rows, characteristic ~ group, value.var = c("pct_unw", "pct_w"))
fwrite(tab, file.path(OUT_DIR, "excluded_vs_included.csv")); print(tab)

## ============================================================ (10) date handling
dh <- pr[f_bp == TRUE, .(n = .N, bio_date_valid = sum(!is.na(bio_date)),
                         bio_date_missing_or_invalid = sum(is.na(bio_date)),
                         fallback_to_interview_date = sum(is.na(bio_date) & !is.na(int_date)))]
fwrite(dh, file.path(OUT_DIR, "date_handling.csv")); print(dh)
reg("DATE_FALLBACK_N", f0(dh$fallback_to_interview_date))
reg("DATE_FALLBACK_PCT", formatC(100 * dh$fallback_to_interview_date / dh$n, format = "f", digits = 2))
log_msg("01 done")
