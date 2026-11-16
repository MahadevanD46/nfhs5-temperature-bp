# =============================================================================
# 04c_additions.R — small additions made at the final check (about 5-10 minutes)
#   1. Numbers of participants per temperature category and per state are counted
#      in the model sample (complete covariates), so they add up to the analytic n.
#   2. Step-by-step denominator of the glucose model (Table 2 footnote).
#   3. Minimally adjusted models: fixed effects only, no covariates (STROBE item 16),
#      for all four outcomes; added to all_estimates.csv (Supplementary Table S10).
# Run after 04_models.R. Safe to run more than once.
# =============================================================================
source("code/B_revision/00_config.R")
a <- load_a()
a <- a[a[, complete.cases(t0, age_grp, hv104, hv270, hv025, shb14a, shb14b, shb14c, sh25, sh26)]]
log_msg("model sample:", f0(nrow(a)))

## ---- 1. counts in the model sample ------------------------------------------
f <- file.path(OUT_DIR, "temperature_categories.csv")
if (file.exists(f)) {
  ct <- fread(f); ct[, c("n", "mean_temp") := NULL]
  ct <- merge(ct, a[, .(n = .N, mean_temp = mean(t0)), by = .(bin = as.character(tbin))], by = "bin")
  fwrite(ct, f); log_msg("category n, total:", f0(ct[outcome == "sbp", sum(n)]))
}
f <- file.path(OUT_DIR, "state_slopes.csv")
if (file.exists(f)) {
  ss <- fread(f); ss[, c("n", "clusters") := NULL]
  ss <- merge(ss, a[, .(n = .N, clusters = uniqueN(hv001)), by = hv024], by = "hv024")
  fwrite(ss, f)
}

## ---- 2. glucose denominator ----------------------------------------------------
n_med   <- a[glu_med %in% 1, .N]
n_medna <- a[is.na(glu_med), .N]
b <- a[glu_med %in% 0]
n_noval <- b[is.na(glu), .N]; b <- b[!is.na(glu)]
n_meal  <- b[is.na(ate_cat), .N]; b <- b[!is.na(ate_cat)]
parts <- c(
  sprintf("users of glucose-lowering medication (%s)", f0(n_med)),
  if (n_medna > 0) sprintf("participants with unknown medication status (%s)", f0(n_medna)),
  sprintf("participants without a valid glucose value (not measured or outside 40–500 mg/dL; %s)", f0(n_noval)),
  if (n_meal > 0) sprintf("participants with missing time since the last meal (%s)", f0(n_meal)))
txt <- paste0("denominator ", f0(nrow(b)), " after excluding, in this order, ",
              paste(parts[-length(parts)], collapse = ", "), " and ", parts[length(parts)])
reg("GLU_DENOM_TEXT", txt); reg("N_GLU_MODEL", f0(nrow(b)))
cat("\n", txt, "\n\n")

## ---- 3. minimally adjusted models (fixed effects only) ---------------------------
a[, htn100 := 100 * htn_measured][, lglu100 := 100 * log(glu)]
LBL0 <- "Fixed effects only, no covariates (minimally adjusted)"
spec <- list(sbp = list("sbp", a), dbp = list("dbp", a), htn = list("htn100", a), glu = list("lglu100", b[, lglu100 := 100 * log(glu)]))
rows <- rbindlist(lapply(names(spec), function(o) {
  m <- feols(as.formula(paste(spec[[o]][[1]], "~ t0 |", FE)), data = spec[[o]][[2]], weights = ~wt,
             cluster = ~hv001, lean = TRUE, mem.clean = TRUE, notes = FALSE)
  r <- row1(LBL0, m, "t0", outcome = o); invisible(gc()); r
}))
print(rows[, .(outcome, n, est, lo, hi)])
for (o in rows$outcome) reg(paste0("MINADJ_", toupper(o), "_CI"), ci(rows[outcome == o, est], rows[outcome == o, lo], rows[outcome == o, hi]))

f <- file.path(OUT_DIR, "all_estimates.csv")
all_est <- fread(f)[model != LBL0]
main_i <- max(which(all_est$model == "Main model"))
all_est <- rbind(all_est[seq_len(main_i)], rows, all_est[-seq_len(main_i)], fill = TRUE)   # directly below the main model
fwrite(all_est, f)
n_tab <- all_est[model == "Main model" & outcome == "glu", n][1]
if (!isTRUE(n_tab == nrow(b))) warning("Glucose model n in all_estimates.csv (", n_tab, ") differs from the step-by-step count (", nrow(b), "): tell the analyst.")
# range of systolic / diastolic estimates across all rows (same rule as in 04_models.R)
for (o in c("sbp", "dbp")) {
  s <- all_est[outcome == o & term %in% c("t0", "t_int", "t_disp", "anom0", "t_mean02") & !grepl("Conley|SEs ", model)]
  reg(paste0("SENS_RANGE_", toupper(o)), paste0(sub("^-", "−", f2(min(s$est))), " to ", sub("^-", "−", f2(max(s$est)))))
}
log_msg("04c done")
