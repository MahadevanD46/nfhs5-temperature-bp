# =============================================================================
# 05_standardisation.R — model-based temperature standardisation and its
#                        sensitivity to the assumptions (comments 21, 22, 36)
#
# For participant i with measured SBP_i, DBP_i, measurement-day temperature T_i
# and reference temperature R_i:
#     SBP*_i = snap( SBP_i - [g_S(T_i) - g_S(R_i)] )
#     DBP*_i = snap( DBP_i - [g_D(T_i) - g_D(R_i)] )
#     elevated*_i = 1{SBP*_i >= 140 or DBP*_i >= 90}
# snap() rounds to the 0.5-mmHg grid. Each participant's own SBP-DBP pair is
# shifted, so the joint distribution (and the SBP-DBP correlation) is kept.
# Main analysis: g linear (g_S(T) = beta_S * T, g_D(T) = beta_D * T), R = annual mean.
# Variants (Supplementary Table S15):
#   UNROUNDED  main without snapping to the 0.5-mmHg grid
#   SEASON     season-specific slopes (by season of measurement)
#   SPLINE     g = restricted cubic spline (script 04)
#   SUPPORT    R clamped to the range of temperatures observed in the district
#   STATE      state-specific empirical-Bayes slopes
#   SEASNORM   R = cluster normal for the season of measurement
#   DOYNORM    R = cluster normal for the calendar day (removes anomalies only)
#   NOALT      clusters <=2,000 m only, slopes re-estimated in that sample
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages(library(Hmisc))
a <- load_a()
a <- a[!is.na(t_annual) & !is.na(t0)]
mv <- readRDS(file.path(OUT_DIR, "models_v2.rds"))
sp <- readRDS(file.path(OUT_DIR, "spline_fit.rds"))
bS <- unname(mv$b["sbp.t0"]); bD <- unname(mv$b["dbp.t0"])
if (is.na(bS)) { bS <- unname(mv$b[[1]]); bD <- unname(mv$b[[2]]) }
reg("BETA_S", formatC(bS, format = "f", digits = 3)); reg("BETA_D", formatC(bD, format = "f", digits = 3))

#' prevalence after shifting by (dS, dD); returns person-level indicator
std_ind <- function(dS, dD, round = TRUE) {
  s <- a$sbp - dS; d <- a$dbp - dD
  if (round) { s <- snap(s); d <- snap(d) }
  as.integer(s >= 140 | d >= 90)
}
lin <- function(b, Tm, R) b * (Tm - R)

rng <- a[, .(tmin = min(t0), tmax = max(t0)), by = shdist]
a <- merge(a, rng, by = "shdist")
R_support <- pmin(pmax(a$t_annual, a$tmin), a$tmax)

sea <- mv$season
gS <- setNames(sea[outcome == "sbp", est], sea[outcome == "sbp", level])
gD <- setNames(sea[outcome == "dbp", est], sea[outcome == "dbp", level])
bS_sea <- gS[as.character(a$season)]; bD_sea <- gD[as.character(a$season)]

stl <- mv$state
bS_st <- stl[outcome == "sbp", setNames(eb, hv024)][as.character(a$hv024)]
bD_st <- stl[outcome == "dbp", setNames(eb, hv024)][as.character(a$hv024)]

kn <- sp$knots
Xs <- function(x) rcspline.eval(x, knots = kn, inclx = TRUE)
dXs <- Xs(a$t0) - Xs(a$t_annual)
spl_S <- as.numeric(dXs %*% sp$coef$sbp$b); spl_D <- as.numeric(dXs %*% sp$coef$dbp$b)
rm(dXs)

V <- list(
  MAIN      = list(lbl = "Main: linear pooled slopes, annual mean",  dS = lin(bS, a$t0, a$t_annual), dD = lin(bD, a$t0, a$t_annual)),
  UNROUNDED = list(lbl = "Main without rounding to 0.5 mmHg",        dS = lin(bS, a$t0, a$t_annual), dD = lin(bD, a$t0, a$t_annual), round = FALSE),
  SEASON    = list(lbl = "Season-specific slopes",                   dS = lin(bS_sea, a$t0, a$t_annual), dD = lin(bD_sea, a$t0, a$t_annual)),
  SPLINE    = list(lbl = "Restricted cubic spline (non-linear)",     dS = spl_S, dD = spl_D),
  SUPPORT   = list(lbl = "Reference restricted to observed district range", dS = lin(bS, a$t0, R_support), dD = lin(bD, a$t0, R_support)),
  STATE     = list(lbl = "State-specific (empirical-Bayes) slopes",  dS = lin(bS_st, a$t0, a$t_annual), dD = lin(bD_st, a$t0, a$t_annual)),
  SEASNORM  = list(lbl = "Reference = season-specific normal",       dS = lin(bS, a$t0, a$t_season_norm), dD = lin(bD, a$t0, a$t_season_norm)),
  DOYNORM   = list(lbl = "Reference = calendar-day normal (anomaly only)", dS = lin(bS, a$t0, a$t_doy), dD = lin(bD, a$t0, a$t_doy)))
for (v in names(V)) a[, (paste0("std_", v)) := std_ind(V[[v]]$dS, V[[v]]$dD, isTRUE(V[[v]]$round %||% TRUE))]

# NOALT: re-estimate slopes in clusters <= 2,000 m
lo <- a[hv040 <= SET$high_alt]
mS <- fit("sbp", "t0", lo); mD <- fit("dbp", "t0", lo)
a[, std_NOALT := NA_integer_]
a[hv040 <= SET$high_alt, std_NOALT := as.integer(
  snap(sbp - coef(mS)["t0"] * (t0 - t_annual)) >= 140 | snap(dbp - coef(mD)["t0"] * (t0 - t_annual)) >= 90)]
V$NOALT <- list(lbl = "Clusters <=2,000 m, slopes re-estimated")
reg("NOALT_BETA_S", formatC(coef(mS)["t0"], format = "f", digits = 3)); rm(lo, mS, mD)

## ---- summaries ------------------------------------------------------------
st_n <- a[, .N, by = hv024][N >= SET$min_state_n, hv024]
long <- rbindlist(lapply(names(V), function(v) {
  col <- paste0("std_", v)
  a[hv024 %in% st_n & !is.na(get(col)), .(variant = v, n = .N,
      obs = 100 * weighted.mean(htn_measured, wt), std = 100 * weighted.mean(get(col), wt)), by = hv024]
}))
long[, change := std - obs]
long[, `:=`(rank_obs = frank(-obs), rank_std = frank(-std)), by = variant]
long[, state := state_lab[as.character(hv024)]]
fwrite(long, file.path(OUT_DIR, "standardisation_variants_long.csv"))
wide <- dcast(long, state ~ variant, value.var = "change")
setcolorder(wide, c("state", names(V)))
fwrite(wide, file.path(OUT_DIR, "standardisation_variants_wide.csv"))

pooled <- rbindlist(lapply(names(V), function(v) {
  col <- paste0("std_", v)
  a[!is.na(get(col)), .(variant = v, label = V[[v]]$lbl, obs = 100 * weighted.mean(htn_measured, wt),
                        std = 100 * weighted.mean(get(col), wt))]
}))
main_ch <- long[variant == "MAIN", .(hv024, main_change = change, main_rank = rank_std)]
summ <- merge(long, main_ch, by = "hv024")[state != "Ladakh", .(
  corr_with_main = cor(change, main_change),
  same_direction = sum(sign(change) == sign(main_change)),
  states = .N, min_change = min(change), max_change = max(change),
  delhi_rank = rank_std[state == "Delhi"], punjab_change = change[state == "Punjab"]), by = variant]
summ <- merge(pooled, summ, by = "variant")
fwrite(summ, file.path(OUT_DIR, "standardisation_variants_summary.csv")); print(summ)
for (v in names(V)) {
  z <- summ[variant == v]
  reg(paste0("STDV_", v, "_POOLED"), f1(z$std))
  reg(paste0("STDV_", v, "_RANGE"), paste0(sub("^-", "\u2212", f1(z$min_change)), " to ", f1(z$max_change)))
  reg(paste0("STDV_", v, "_COR"), formatC(z$corr_with_main, format = "f", digits = 2))
  reg(paste0("STDV_", v, "_DELHI"), z$delhi_rank)
}
reg("STDV_POOLED_RANGE", paste0(f1(min(summ[variant != "UNROUNDED", std])), "% to ",
                                 f1(max(summ[variant != "UNROUNDED", std])), "%"))
reg("STDV_COR_MIN", formatC(min(summ[variant != "MAIN", corr_with_main]), format = "f", digits = 2))

# district-level observed and standardised prevalence for LISA (script 07)
dist <- a[, .(n = .N, hv024 = hv024[1], obs = 100 * weighted.mean(htn_measured, wt),
              std = 100 * weighted.mean(std_MAIN, wt)), by = shdist]
fwrite(dist, file.path(OUT_DIR, "district_prevalence.csv"))
log_msg("05 done")
