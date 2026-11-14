# =============================================================================
# 00_config.R — shared settings and helpers for the REVISION analyses
#
# Run order:  code/A_original_pipeline/nfhs5_temperature_full_period.R
#             then code/B_revision/run_all.R  (or scripts 01-08 one by one)
# Results are printed to the console and written as plain CSV files to output_revision/.
# Working directory = your project folder, i.e. the folder that contains
#   IAPR7EFL.DTA, IAGE7AFL.shp, climate/ and output_full/
# Every script starts with  source("code/B_revision/00_config.R")
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(fixest) })
# Make sure R runs in UTF-8 (default on Windows from R 4.2 and on macOS).
if (!isTRUE(l10n_info()[["UTF-8"]]))
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", if (.Platform$OS.type == "windows") ".UTF-8" else "C.UTF-8")))

## ----------------------------------------------------------------- EDIT HERE
PR_FILE   <- "IAPR7EFL.DTA"
GE_FILE   <- "IAGE7AFL.shp"
CLIM_DIR  <- "climate"
IN_DIR    <- "output_full"          # outputs of the original pipeline
OUT_DIR   <- "output_revision"      # everything new goes here
FIG_DIR   <- file.path(OUT_DIR, "figures")
# DHS Spatial Data Repository district boundaries (NFHS-5, admin2) and the
# attribute that equals shdist. Open the .shp in R (sf::st_read) to check.
DIST_SHP      <- "boundaries/sdr_subnational_boundaries2.shp"
DIST_ID_FIELD <- "REGCODE"
# Optional IMD 1-degree gridded daily Tmax/Tmin for independent validation
IMD_DIR   <- "imd"
## ---------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE); dir.create(FIG_DIR, showWarnings = FALSE)

SET <- list(
  seed = 20261112,
  n_boot = 500,                 # bootstrap replicates; 200 is acceptable for a first run
  conley_km = c(50, 100),
  spline_pct = c(.05, .35, .65, .95),
  spline_ref = 22.5,            # centre of the 20-25 degC reference category
  min_state_n = 500, min_district_n = 100,
  lisa_k = 6, lisa_nsim = 999, lisa_alpha = 0.05,
  high_alt = 2000,
  n_disp_draws = 200
)
set.seed(SET$seed)
setFixest_nthreads(max(1, parallel::detectCores() - 1))

state_lab <- c("1"="Jammu & Kashmir","2"="Himachal Pradesh","3"="Punjab","4"="Chandigarh",
  "5"="Uttarakhand","6"="Haryana","7"="Delhi","8"="Rajasthan","9"="Uttar Pradesh","10"="Bihar",
  "11"="Sikkim","12"="Arunachal Pradesh","13"="Nagaland","14"="Manipur","15"="Mizoram",
  "16"="Tripura","17"="Meghalaya","18"="Assam","19"="West Bengal","20"="Jharkhand","21"="Odisha",
  "22"="Chhattisgarh","23"="Madhya Pradesh","24"="Gujarat","25"="Dadra & Nagar Haveli and Daman & Diu",
  "27"="Maharashtra","28"="Andhra Pradesh","29"="Karnataka","30"="Goa","31"="Lakshadweep",
  "32"="Kerala","33"="Tamil Nadu","34"="Puducherry","35"="Andaman & Nicobar Islands",
  "36"="Telangana","37"="Ladakh")
# NFHS regions (used for climatic-zone heterogeneity, comment 36)
region_of <- c("1"="North","2"="North","3"="North","4"="North","5"="North","6"="North",
  "7"="North","8"="North","37"="North","9"="Central","22"="Central","23"="Central",
  "10"="East","19"="East","20"="East","21"="East","11"="Northeast","12"="Northeast",
  "13"="Northeast","14"="Northeast","15"="Northeast","16"="Northeast","17"="Northeast",
  "18"="Northeast","24"="West","25"="West","27"="West","30"="West","28"="South",
  "29"="South","31"="South","32"="South","33"="South","34"="South","35"="South","36"="South")

## ----------------------------------------------------- model specification
## Identical to the original pipeline (Section 7 of nfhs5_temperature_full_period.R)
CTRL  <- "age_grp*hv104 + i(hv270) + hv025 + i(shb14a) + i(shb14b) + i(shb14c) + i(sh25) + i(sh26)"
FE    <- "shdist + ym + dow"
FE_DM <- "shdist^ym + dow"

#' fixest formula  lhs ~ x + CTRL | fe
F <- function(lhs, x = "t0", fe = FE, extra = NULL)
  as.formula(paste(lhs, "~", paste(c(x, extra, CTRL), collapse = " + "), "|", fe))

#' Survey-weighted feols with cluster-robust SEs (PSU = hv001), memory-lean
fit <- function(lhs, x = "t0", data, fe = FE, extra = NULL, cl = ~hv001, sub = NULL,
                vcov = NULL, lean = TRUE) {
  args <- list(fml = F(lhs, x, fe, extra), data = data, weights = ~wt,
               lean = lean, mem.clean = TRUE, notes = FALSE)
  if (is.null(vcov)) args$cluster <- cl else args$vcov <- vcov
  if (!is.null(sub)) args$subset <- sub
  m <- do.call(feols, args); invisible(gc()); m
}

#' One coefficient as a tidy row (scale = 100 for pp or %)
row1 <- function(label, m, term = "t0", scale = 1, outcome = NA) {
  ct <- coeftable(m)
  idx <- if (term %in% rownames(ct)) term else grep(term, rownames(ct), fixed = TRUE)[1]
  b <- ct[idx, 1]; se <- ct[idx, 2]
  data.table(outcome = outcome, model = label, term = if (is.character(idx)) idx else rownames(ct)[idx], n = nobs(m),
             est = scale * b, lo = scale * (b - 1.959964 * se), hi = scale * (b + 1.959964 * se),
             se = scale * se, p = 2 * pnorm(-abs(b / se)))
}

#' Wald test that a set of coefficients are all zero (regex on names)
wald_p <- function(m, pattern) {
  b <- coef(m); V <- vcov(m); k <- grep(pattern, names(b))
  W <- as.numeric(t(b[k]) %*% solve(V[k, k]) %*% b[k])
  c(chi2 = W, df = length(k), p = pchisq(W, length(k), lower.tail = FALSE))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
snap <- function(x) round(2 * x) / 2               # 0.5-mmHg measurement grid
wmean <- function(x, w) { ok <- !is.na(x); sum(x[ok] * w[ok]) / sum(w[ok]) }
wcor  <- function(x, y, w) { ok <- complete.cases(x, y, w)
  stats::cov.wt(cbind(x[ok], y[ok]), wt = w[ok], cor = TRUE)$cor[1, 2] }

## ----------------------------------------------------- results to console
#' Print a labelled result to the console
reg <- function(key, value, note = "") { cat(sprintf("  %-28s %s\n", key, value)); invisible(value) }
f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
f2 <- function(x) formatC(x, format = "f", digits = 2)
pct <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
fp <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
pk <- function(p) ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, format = "f", digits = 3)))
ci <- function(est, lo, hi, d = 2) {
  f <- function(x) sub("^-", "\u2212", formatC(x, format = "f", digits = d))
  paste0(f(est), " (", f(lo), " to ", f(hi), ")")
}
#' print estimate and 95% CI of a tidy row
reg_row <- function(key, r, d = 2) reg(key, ci(r$est, r$lo, r$hi, d))
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

#' Analysis data from the original pipeline + new variables (script 02)
load_a <- function(extra = TRUE, cols = NULL) {
  f <- file.path(OUT_DIR, "analysis_dataset_v2.rds")
  if (!extra || !file.exists(f)) f <- file.path(IN_DIR, "analysis_dataset.rds")
  a <- readRDS(f); setDT(a)
  if (!is.null(cols)) a <- a[, intersect(cols, names(a)), with = FALSE]
  a
}
cat("Revision config loaded | fixest", as.character(packageVersion("fixest")), "\n")
