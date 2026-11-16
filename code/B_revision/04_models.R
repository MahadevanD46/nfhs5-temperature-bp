# =============================================================================
# 04_models.R — all regression analyses requested in the review
#   main models + derived contrasts with CIs               (comments 4, 25)
#   sensitivity analyses for ALL outcomes                  (9, 10, 12, 15, 16, 34, 35)
#   alternative variance estimators incl. Conley           (16)
#   temperature categories with n                          (41)
#   restricted cubic splines + test of non-linearity       (17)
#   pooled interaction models + Wald tests                 (18, 34, 36)
#   placebo diagnostics for all outcomes                   (26)
#   state- and region-specific slopes, heterogeneity       (36)
# Outputs: output_revision/*.csv (one file per supplementary table) and
#          output_revision/models_v2.rds (coefficients and vcov only)
# Memory: models are fitted lean and discarded; ~8-16 GB RAM is sufficient.
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages({ library(Hmisc); library(metafor); library(spdep) })
a <- load_a()
# model sample (complete covariates), so that all counts add up to the analytic n
a <- a[a[, complete.cases(t0, age_grp, hv104, hv270, hv025, shb14a, shb14b, shb14c, sh25, sh26)]]
a[, htn100 := 100 * htn_measured][, htnrep100 := 100 * htn_report]
a[, lglu100 := 100 * log(glu)]                       # coefficient = % change (approx.)

OUT <- list(
  sbp = list(lhs = "sbp",       x0 = NULL,          sub = NULL),
  dbp = list(lhs = "dbp",       x0 = NULL,          sub = NULL),
  htn = list(lhs = "htn100",    x0 = NULL,          sub = NULL),
  glu = list(lhs = "lglu100",   x0 = "i(ate_cat)",  sub = quote(glu_med == 0)))
LBL <- c(sbp = "Systolic BP (mmHg)", dbp = "Diastolic BP (mmHg)",
         htn = "Measured elevated BP (pp)", glu = "Random glucose (%)")

#' Fit outcome o with exposure x; extra terms, subset expression, FE, clustering
run <- function(o, x = "t0", extra = NULL, sub = NULL, fe = FE, cl = ~hv001, vcov = NULL,
                lhs = NULL) {
  sp <- OUT[[o]]
  s <- sp$sub
  if (!is.null(sub)) s <- if (is.null(s)) sub else bquote(.(s) & .(sub))
  d <- if (is.null(s)) a else a[eval(s, a)]
  fit(if (is.null(lhs)) sp$lhs else lhs, x, d, fe = fe, extra = c(sp$x0, extra), cl = cl, vcov = vcov)
}
res <- list()
keep <- function(lbl, o, m, term = "t0") {
  r <- row1(lbl, m, term, outcome = o); res[[length(res) + 1]] <<- r; r }

## =================================================== main + contrasts (4, 25)
main <- list()
for (o in names(OUT)) { main[[o]] <- run(o); keep("Main model", o, main[[o]]) }
p5 <- quantile(a$t0, .05); p95 <- quantile(a$t0, .95); q1 <- quantile(a$t0, .25); q3 <- quantile(a$t0, .75)
reg("T_P5", f1(p5)); reg("T_P95", f1(p95)); reg("T_W595", f1(p95 - p5)); reg("T_WIQR", f1(q3 - q1))
for (o in c("sbp", "dbp", "htn")) {
  r <- row1("", main[[o]])
  for (w in list(c("C595", p95 - p5), c("CIQR", q3 - q1))) {
    k <- as.numeric(w[2])
    reg(paste0(w[1], "_", toupper(o)), ci(-r$est * k, -r$hi * k, -r$lo * k, 1))  # cold minus warm
  }
}

## =================================================== sensitivity analyses
SENS <- list(
  list(id = "INTDATE",  lbl = "Exposure on interview date",            x = "t_int"),
  list(id = "SAMEDAY",  lbl = "Measured on interview day only",        sub = quote(same_day == 1)),
  list(id = "GAP1",     lbl = "Interview-biomarker gap <=1 day",       sub = quote(abs(gap) <= 1)),
  list(id = "NOFALLBK", lbl = "Excluding interview-date fallback",     sub = quote(date_fallback == 0)),
  list(id = "DISP",     lbl = "Displacement-weighted temperature",     x = "t_disp"),
  list(id = "HUMID",    lbl = "Adjusted for humidity and precipitation", extra = c("rh", "i(tp_cat)")),
  list(id = "ANOM",     lbl = "Temperature anomaly",                   x = "anom0"),
  list(id = "MEAN02",   lbl = "Mean temperature, days 0-2",            x = "t_mean02"),
  list(id = "NOMEDS",   lbl = "Excluding antihypertensive users",      sub = quote(bp_med == 0)),
  list(id = "NOALT",    lbl = "Excluding clusters >2,000 m",           sub = quote(hv040 <= 2000)),
  list(id = "NOFILL",   lbl = "Excluding coastal-infilled clusters",   sub = quote(filled == FALSE)),
  list(id = "NOALTFILL",lbl = "Excluding >2,000 m and coastal-infilled", sub = quote(hv040 <= 2000 & filled == FALSE)),
  list(id = "BMI",      lbl = "Adjusted for BMI (BMI subsample)",      extra = "bmi", sub = quote(sub_bmi == TRUE)),
  list(id = "DM",       lbl = "District-by-month fixed effects",       fe = FE_DM),
  list(id = "HOUR",     lbl = "Hour-of-measurement fixed effects",     extra = "i(hour_cat)"),
  list(id = "INTV",     lbl = "Interviewer (within state) fixed effects", fe = paste(FE, "+ hv024^hv018")),
  list(id = "CLDIST",   lbl = "SEs clustered by district",             cl = ~shdist),
  list(id = "CL2WAY",   lbl = "SEs two-way clustered (cluster, date)", cl = ~hv001 + date),
  list(id = "CLCELL",   lbl = "SEs clustered by ERA5-Land cell x date", cl = ~cell_date),
  list(id = "TRIM",     lbl = "Excluding top/bottom 1% of residual temperature", sub = quote(trim_ok == TRUE)))

r1 <- feols(t0 ~ 1 | shdist + ym + dow, a, weights = ~wt, notes = FALSE)
a[, t_res := resid(r1, na.rm = FALSE)]
lim <- quantile(a$t_res, c(.01, .99), na.rm = TRUE)
a[, trim_ok := !is.na(t_res) & t_res >= lim[1] & t_res <= lim[2]]
rm(r1); invisible(gc())

for (s in SENS) for (o in names(OUT)) {
  if (s$id == "NOMEDS" && o == "glu") next
  if (s$id == "HUMID" && !all(c("rh", "tp_cat") %in% names(a))) next
  m <- tryCatch(run(o, x = if (is.null(s$x)) "t0" else s$x, extra = s$extra, sub = s$sub,
                    fe = if (is.null(s$fe)) FE else s$fe, cl = if (is.null(s$cl)) ~hv001 else s$cl),
                error = function(e) { log_msg(s$id, o, "FAILED:", conditionMessage(e)); NULL })
  if (is.null(m)) next
  r <- keep(s$lbl, o, m, term = if (is.null(s$x)) "t0" else s$x)
  reg_row(paste(s$id, toupper(o), sep = "_"), r)
  if (o == "sbp") reg(paste0(s$id, "_N"), f0(r$n))
  rm(m); invisible(gc())
}
# composite outcome (comment 9): elevated BP OR medication
m <- run("htn", lhs = "htnrep100"); r <- keep("Composite: >=140/90 mmHg or medication", "htn", m)
reg_row("COMPOSITE", r); rm(m)

# normalised weights give identical estimates and SEs (comment 16)
a[, wt_norm := wt / mean(wt)]
mw <- feols(F("sbp"), a, weights = ~wt_norm, cluster = ~hv001, lean = TRUE, notes = FALSE)
reg("WT_NORM_CHECK", sprintf("%.6f (SE %.6f) vs %.6f (SE %.6f)", coef(mw)["t0"], se(mw)["t0"],
                             coef(main$sbp)["t0"], se(main$sbp)["t0"]))
rm(mw)

## =================================================== Conley SEs (16)
# Spatial HAC with a uniform kernel, computed from cluster-aggregated scores
# via the Frisch-Waugh-Lovell theorem (exact for the t0 coefficient).
conley <- function(o, cutoffs) {
  sp <- OUT[[o]]; d <- if (is.null(sp$sub)) a else a[eval(sp$sub, a)]
  m <- feols(F(sp$lhs, extra = sp$x0), d, weights = ~wt, notes = FALSE, lean = FALSE)
  e <- resid(m, na.rm = FALSE); b <- coef(m)["t0"]
  rx <- feols(F("t0", x = "1", extra = sp$x0), d[obs(m)], weights = ~wt, notes = FALSE)
  # NB: F("t0", x = "1") regresses t0 on the covariates and FE
  xt <- resid(rx)
  dd <- d[obs(m), .(hv001, lat, lon, w = wt)][, `:=`(s = w * xt * e[obs(m)])]
  D <- sum(dd$w * xt^2)
  S <- dd[, .(S = sum(s), lat = lat[1], lon = lon[1]), by = hv001][!is.na(lat)]
  out <- data.table(cutoff_km = 0, se = sqrt(sum(S$S^2)) / D)
  for (k in cutoffs) {
    nb <- dnearneigh(cbind(S$lon, S$lat), 0, k, longlat = TRUE)
    cross <- sum(S$S * vapply(nb, function(j) if (length(j) && j[1] > 0) sum(S$S[j]) else 0, 0))
    out <- rbind(out, data.table(cutoff_km = k, se = sqrt(sum(S$S^2) + cross) / D))
  }
  out[, `:=`(outcome = o, est = b, beta_fwl = sum(dd$w * xt * d[obs(m), get(sp$lhs)]) / D,
             se_fixest_cluster = se(summary(m, cluster = ~hv001))["t0"])]
  rm(m, rx); invisible(gc()); out
}
cz <- rbindlist(lapply(names(OUT), conley, cutoffs = SET$conley_km))
fwrite(cz, file.path(OUT_DIR, "conley_se.csv")); print(cz)
for (o in names(OUT)) for (k in SET$conley_km) {
  z <- cz[outcome == o & cutoff_km == k]
  res[[length(res) + 1]] <- data.table(outcome = o, model = sprintf("Conley SEs, %d km", k), term = "t0",
    n = nobs(main[[o]]), est = z$est * 1, lo = z$est - 1.959964 * z$se, hi = z$est + 1.959964 * z$se,
    se = z$se, p = 2 * pnorm(-abs(z$est / z$se)))
  reg(sprintf("CONLEY%d_%s_CI", k, toupper(o)), ci(z$est, z$est - 1.96 * z$se, z$est + 1.96 * z$se))
}

## =================================================== categories (41)
cat_tab <- rbindlist(lapply(c("sbp", "dbp", "htn"), function(o) {
  m <- run(o, x = "i(tbin, ref = '20-25')")
  ct <- as.data.table(coeftable(m), keep.rownames = "term")[grepl("^tbin::", term)]
  ct[, bin := sub("tbin::", "", term)]
  out <- rbind(ct[, .(bin, est = Estimate, lo = Estimate - 1.959964 * `Std. Error`,
                      hi = Estimate + 1.959964 * `Std. Error`)],
               data.table(bin = "20-25", est = 0, lo = NA, hi = NA))
  out[, outcome := o]
}))
cat_n <- a[, .(n = .N, mean_temp = mean(t0)), by = .(bin = as.character(tbin))]
cat_tab <- merge(cat_tab, cat_n, by = "bin")
fwrite(cat_tab, file.path(OUT_DIR, "temperature_categories.csv"))

## =================================================== splines (17)
kn <- quantile(a$t0, SET$spline_pct)
X <- rcspline.eval(a$t0, knots = kn, inclx = TRUE)
a[, c("t_rcs1", "t_rcs2", "t_rcs3") := as.data.table(X)]
grid <- seq(quantile(a$t0, .01), quantile(a$t0, .99), by = 0.25)
Xg <- rcspline.eval(grid, knots = kn, inclx = TRUE); Xr <- rcspline.eval(SET$spline_ref, knots = kn, inclx = TRUE)
spl_curves <- list(); spl_tests <- list(); spl_coef <- list()
for (o in names(OUT)) {
  m <- run(o, x = "t_rcs1 + t_rcs2 + t_rcs3")
  b <- coef(m)[c("t_rcs1", "t_rcs2", "t_rcs3")]; V <- vcov(m)[names(b), names(b)]
  dX <- sweep(Xg, 2, Xr)
  fit_ <- dX %*% b; se_ <- sqrt(rowSums((dX %*% V) * dX))
  spl_curves[[o]] <- data.table(outcome = o, temp = grid, est = fit_[, 1],
                                lo = fit_[, 1] - 1.96 * se_, hi = fit_[, 1] + 1.96 * se_)
  w <- wald_p(m, "^t_rcs[23]$")
  spl_tests[[o]] <- data.table(outcome = o, chi2 = w["chi2"], df = w["df"], p_nonlinear = w["p"])
  spl_coef[[o]] <- list(b = b, V = V)                   # used by 05_standardisation.R
  reg(paste0("SPLINE_P_", toupper(o)), pk(w["p"]))
  reg(paste0("SPLINE_CHI_", toupper(o)), f1(w["chi2"]))
  rm(m); invisible(gc())
}
fwrite(rbindlist(spl_curves), file.path(OUT_DIR, "spline_curves.csv"))
fwrite(rbindlist(spl_tests), file.path(OUT_DIR, "spline_tests.csv"))
saveRDS(list(knots = kn, coef = spl_coef), file.path(OUT_DIR, "spline_fit.rds"))
reg("SPLINE_KNOTS", paste(f1(kn), collapse = ", "))

## =================================================== effect modification (18, 34, 36)
# fixest cannot parse a ref level containing commas/brackets, e.g. "(69,Inf]", so
# modifiers with such labels get safe integer codes; labels are restored below.
a[, `:=`(age_k = as.integer(age_grp), hour_k = as.integer(factor(hour_cat)))]
LAB <- list(age_k = setNames(levels(a$age_grp), seq_along(levels(a$age_grp))),
            hour_k = setNames(levels(factor(a$hour_cat)), seq_along(levels(factor(a$hour_cat)))))
EM <- list(
  sex    = list(v = "hv104",    lbl = "Sex",      ref = 1),
  age    = list(v = "age_k",    lbl = "Age group", ref = NULL),
  season = list(v = "season",   lbl = "Season",   ref = "winter"),
  hour   = list(v = "hour_k",   lbl = "Hour of measurement", ref = NULL, extra = "i(hour_k)"),
  region = list(v = "region",   lbl = "Region",   ref = "North"),
  phase  = list(v = "shphase",  lbl = "Survey phase", ref = NULL))
em_slopes <- list(); em_tests <- list()
for (nm in names(EM)) for (o in names(OUT)) {
  e <- EM[[nm]]
  ref <- if (is.null(e$ref)) names(sort(-table(a[[e$v]])))[1] else e$ref
  # (i) separate slopes per level, from one pooled model
  m1 <- run(o, x = sprintf("i(%s, t0)", e$v), extra = e$extra)
  ct <- as.data.table(coeftable(m1), keep.rownames = "term")[grepl(paste0("^", e$v, "::.*:t0$"), term)]
  em_slopes[[paste(nm, o)]] <- ct[, .(modifier = e$lbl, outcome = o,
    level = sub(paste0("^", e$v, "::(.*):t0$"), "\\1", term), est = Estimate,
    lo = Estimate - 1.959964 * `Std. Error`, hi = Estimate + 1.959964 * `Std. Error`)]
  # (ii) common slope + interactions -> joint Wald test of equal slopes
  m2 <- run(o, x = sprintf("t0 + i(%s, t0, ref = '%s')", e$v, ref), extra = e$extra)
  # (fixest names a single interaction column "i(factor_var = v, ...)" instead of "v::level:t0")
  w <- wald_p(m2, paste0("^", e$v, "::.*:t0$|factor_var = ", e$v, ","))   # interactions only
  em_tests[[paste(nm, o)]] <- data.table(modifier = e$lbl, outcome = o, reference = ref,
                                         chi2 = w["chi2"], df = w["df"], p = w["p"])
  reg(sprintf("EM_%s_P_%s", toupper(nm), toupper(o)), pk(w["p"]))
  reg(sprintf("EM_%s_CHI_%s", toupper(nm), toupper(o)), paste0(f1(w["chi2"]), ", ", w["df"], " df"))
  rm(m1, m2); invisible(gc())
}
em_s <- rbindlist(em_slopes); em_t <- rbindlist(em_tests)
for (v in names(LAB)) {
  lbl <- EM[[which(sapply(EM, `[[`, "v") == v)]]$lbl
  em_s[modifier == lbl, level := LAB[[v]][level]]
  em_t[modifier == lbl, reference := LAB[[v]][reference]]
}
fwrite(em_s, file.path(OUT_DIR, "effect_modification_slopes.csv"))
fwrite(em_t, file.path(OUT_DIR, "effect_modification_tests.csv"))
age_sbp <- em_s[modifier == "Age group" & outcome == "sbp"]
reg("AGE_SLOPES_SBP", paste(sprintf("%s: %s", age_sbp$level, ci(age_sbp$est, age_sbp$lo, age_sbp$hi)), collapse = "; "))
# fan / cooler (as in the original analysis)
for (o in names(OUT)) {
  m <- run(o, x = "t0 * fan + t0 * cooler")
  keep("Interaction: temperature x fan", o, m, "t0:fan"); keep("Interaction: temperature x cooler", o, m, "t0:cooler")
}
# season slopes after adjustment for humidity and rain (is the monsoon slope humidity?)
if (all(c("rh", "tp_cat") %in% names(a))) {
  m <- run("sbp", x = "i(season, t0)", extra = c("rh", "i(tp_cat)"))
  ms <- as.data.table(coeftable(m), keep.rownames = "term")[grepl("^season::", term)]
  ms[, level := sub("^season::(.*):t0$", "\\1", term)]
  fwrite(ms, file.path(OUT_DIR, "season_slopes_humidity_adjusted.csv"))
  z <- ms[level == "monsoon"]
  reg("MONSOON_HUMID_CI", ci(z$Estimate, z$Estimate - 1.96 * z$`Std. Error`, z$Estimate + 1.96 * z$`Std. Error`))
  m2 <- run("sbp", x = "t0 + i(season, t0, ref = 'winter')", extra = c("rh", "i(tp_cat)"))
  reg("SEASON_HUMID_P", pk(wald_p(m2, "^season::.*:t0$")["p"])); rm(m, m2)
}

## =================================================== placebo diagnostics (26)
rx <- feols(c(t0, t_lead7) ~ 1 | shdist + ym + dow, a, weights = ~wt, notes = FALSE)
rho_raw <- wcor(a$t0, a$t_lead7, a$wt)
rho_res <- wcor(resid(rx[[1]], na.rm = FALSE), resid(rx[[2]], na.rm = FALSE), a$wt)
reg("PLACEBO_COR", formatC(rho_raw, format = "f", digits = 2))
reg("PLACEBO_RESCOR", formatC(rho_res, format = "f", digits = 2))
reg("PLACEBO_VIF", f2(1 / (1 - rho_res^2)))
pl <- rbindlist(lapply(names(OUT), function(o) {
  ml <- run(o, x = "t_lead7"); mj <- run(o, x = "t0 + t_lead7")
  out <- rbind(row1("Lead alone", ml, "t_lead7", outcome = o),
               row1("Joint: same day", mj, "t0", outcome = o),
               row1("Joint: lead", mj, "t_lead7", outcome = o))
  out[, se_ratio_joint_vs_main := c(NA, se(mj)["t0"] / se(main[[o]])["t0"], NA)]
  for (k in 1:3) reg_row(sprintf("PL%d_%s", k, toupper(o)), out[k])
  out
}))
fwrite(pl, file.path(OUT_DIR, "placebo_all_outcomes.csv"))
reg("PLACEBO_SEINF", f2(pl[outcome == "sbp" & model == "Joint: same day", se_ratio_joint_vs_main]))
rm(rx); invisible(gc())

## =================================================== heterogeneity by state (36)
st_slopes <- list()
for (o in c("sbp", "dbp")) {
  m <- run(o, x = "i(hv024, t0)")
  ct <- as.data.table(coeftable(m), keep.rownames = "term")[grepl("^hv024::", term)]
  ct[, hv024 := as.integer(sub("^hv024::(\\d+):t0$", "\\1", term))]
  st <- ct[, .(outcome = o, hv024, state = state_lab[as.character(hv024)],
               est = Estimate, se = `Std. Error`)]
  st <- merge(st, a[, .(n = .N, clusters = uniqueN(hv001)), by = hv024], by = "hv024")
  st <- st[n >= SET$min_state_n & clusters >= 20]       # as in Table 3; avoids unstable slopes
  ma <- rma(yi = est, sei = se, data = st, method = "REML")
  st[, eb := as.numeric(blup(ma)$pred)]                 # empirical-Bayes slopes for 05
  st_slopes[[o]] <- st
  reg(paste0("STATE_I2_", toupper(o)), f1(ma$I2)); reg(paste0("STATE_TAU_", toupper(o)), f2(sqrt(ma$tau2)))
  reg(paste0("STATE_QP_", toupper(o)), pk(ma$QEp)); reg(paste0("STATE_POOLED_", toupper(o)),
      ci(ma$b[1], ma$ci.lb, ma$ci.ub))
  reg(paste0("STATE_RANGE_", toupper(o)), paste0(f2(min(st$est)), " to ", f2(max(st$est))))
  reg(paste0("STATE_NEG_", toupper(o)), st[est < 0, .N])
  reg(paste0("STATE_PI_", toupper(o)), {p <- predict(ma); paste0(f2(p$pi.lb), " to ", f2(p$pi.ub))})
  rm(m); invisible(gc())
}
fwrite(rbindlist(st_slopes), file.path(OUT_DIR, "state_slopes.csv"))

# leave-one-state-out (is the estimate driven by a few places?)
loso <- rbindlist(lapply(sort(unique(a$hv024)), function(s) {
  m <- fit("sbp", "t0", a[hv024 != s])
  data.table(state_dropped = state_lab[as.character(s)], est = coef(m)["t0"], se = se(m)["t0"])
}))
fwrite(loso, file.path(OUT_DIR, "leave_one_state_out.csv"))
reg("LOSO_RANGE", paste0(f2(min(loso$est)), " to ", f2(max(loso$est))))

## =================================================== save
all_est <- rbindlist(res, fill = TRUE)
se_s <- all_est[outcome == "sbp" & term %in% c("t0", "t_int", "t_disp", "anom0", "t_mean02") &
                !grepl("Conley|SEs ", model)]
reg("SENS_RANGE_SBP", paste0(sub("^-", "\u2212", f2(min(se_s$est))), " to ", sub("^-", "\u2212", f2(max(se_s$est)))))
se_d <- all_est[outcome == "dbp" & term %in% c("t0", "t_int", "t_disp", "anom0", "t_mean02") & !grepl("Conley|SEs ", model)]
reg("SENS_RANGE_DBP", paste0(sub("^-", "\u2212", f2(min(se_d$est))), " to ", sub("^-", "\u2212", f2(max(se_d$est)))))
fwrite(all_est, file.path(OUT_DIR, "all_estimates.csv"))
saveRDS(list(b = sapply(main, function(m) coef(m)["t0"]), V = lapply(main, vcov),
             season = em_s[modifier == "Season"], state = rbindlist(st_slopes)),
        file.path(OUT_DIR, "models_v2.rds"))
log_msg("04 done")
