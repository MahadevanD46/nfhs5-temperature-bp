# =============================================================================
# 06_bootstrap.R — design-based uncertainty for standardised prevalence, changes,
#                  ranks and district hotspots (comments 25, 27, 28, 39)
#
# Rao-Wu rescaling bootstrap (Rao, Wu & Yue 1992) for stratified multistage
# designs: in each replicate, n_h - 1 PSUs are drawn with replacement within
# each stratum h (hv022) and weights are multiplied by m_hi * n_h / (n_h - 1).
# Strata with a single PSU are collapsed within district. In EVERY replicate
# both slopes (SBP and DBP) are re-estimated with the full fixed-effects model,
# the standardisation is redone, and state/district prevalence and ranks are
# recomputed — so the intervals combine sampling uncertainty in prevalence AND
# uncertainty in the slopes, including their correlation.
#
# Runtime: about 10-20 s per replicate on 1.6 M rows (8 cores) -> 500 replicates
# take roughly 2-3 hours. The script is resumable: re-run it and it continues
# from the last saved replicate.
# Outputs: output_revision/boot_reps.rds, table3_bootstrap.csv, boot_pooled.csv
# =============================================================================
source("code/B_revision/00_config.R")
a <- load_a(cols = c("hv001", "hv021", "hv022", "hv024", "shdist", "ym", "dow", "wt", "sbp", "dbp",
  "htn_measured", "t0", "t_annual", "age_grp", "hv104", "hv270", "hv025", "shb14a",
  "shb14b", "shb14c", "sh25", "sh26"))
a <- a[!is.na(t_annual) & complete.cases(a)]
invisible(gc())

## ---- design: PSUs within (collapsed) strata ------------------------------------
psu <- unique(a[, .(hv021, hv022, shdist)])
psu[, n_h0 := .N, by = hv022]
psu[, stratum := paste0("s", hv022)]
# strata with a single PSU: pool all such PSUs of a district into one stratum ...
psu[n_h0 == 1, stratum := paste0("d", shdist)]
# ... and if that pooled stratum still has one PSU, add it to the district's largest stratum
psu[, n_tmp := .N, by = stratum]
big <- psu[n_h0 > 1, .N, by = .(shdist, hv022)][order(-N)][, .SD[1], by = shdist]
psu[n_tmp == 1, stratum := fcoalesce(paste0("s", big$hv022[match(shdist, big$shdist)]), "lone_all")]
psu[stratum == "sNA", stratum := "lone_all"]
psu[, n_h := .N, by = stratum]
reg("BOOT_STRATA", f0(uniqueN(psu$stratum))); reg("BOOT_PSU", f0(nrow(psu)))
a[, psu_row := match(hv021, psu$hv021)]

st_keep <- a[, .N, by = hv024][N >= SET$min_state_n, sort(hv024)]
dist_keep <- a[hv024 != 37, .N, by = shdist][N >= SET$min_district_n, sort(shdist)]
fml <- as.formula(paste("c(sbp, dbp) ~ t0 +", CTRL, "|", FE))

one_rep <- function(factor_psu) {
  a[, wb := wt * factor_psu[psu_row]]
  d <- a[wb > 0]
  m <- feols(fml, d, weights = ~wb, lean = TRUE, notes = FALSE, mem.clean = TRUE)
  bS <- coef(m[[1]])["t0"]; bD <- coef(m[[2]])["t0"]
  sh <- d$t0 - d$t_annual
  d[, h_std := as.integer(snap(sbp - bS * sh) >= 140 | snap(dbp - bD * sh) >= 90)]
  st <- d[hv024 %in% st_keep, .(obs = 100 * weighted.mean(htn_measured, wb),
                                std = 100 * weighted.mean(h_std, wb)), keyby = hv024]
  ds <- d[shdist %in% dist_keep, .(obs = 100 * weighted.mean(htn_measured, wb),
                                   std = 100 * weighted.mean(h_std, wb)), keyby = shdist]
  list(bS = bS, bD = bD,
       pooled = d[, c(obs = 100 * weighted.mean(htn_measured, wb), std = 100 * weighted.mean(h_std, wb))],
       st_obs = setNames(st$obs, st$hv024), st_std = setNames(st$std, st$hv024),
       d_obs = setNames(ds$obs, ds$shdist), d_std = setNames(ds$std, ds$shdist))
}

## ---- point estimate (factor = 1) and replicates ---------------------------------
bfile <- file.path(OUT_DIR, "boot_reps.rds")
B <- if (file.exists(bfile)) readRDS(bfile) else list(point = one_rep(rep(1, nrow(psu))), reps = list(),
                                                     rng = NULL)
if (!is.null(B$rng)) .Random.seed <- B$rng
while (length(B$reps) < SET$n_boot) {
  f <- numeric(nrow(psu))
  for (h in split(seq_len(nrow(psu)), psu$stratum)) {
    nh <- length(h)
    m_hi <- tabulate(sample.int(nh, nh - 1, replace = TRUE), nh)
    f[h] <- m_hi * nh / (nh - 1)
  }
  B$reps[[length(B$reps) + 1]] <- one_rep(f)
  B$rng <- .Random.seed
  if (length(B$reps) %% 10 == 0) { saveRDS(B, bfile); log_msg("replicate", length(B$reps)) }
}
saveRDS(B, bfile)

## ---- summaries -------------------------------------------------------------------
R <- B$reps; P <- B$point
qq <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
mat <- function(nm) do.call(cbind, lapply(R, `[[`, nm))
obsM <- mat("st_obs"); stdM <- mat("st_std"); chM <- stdM - obsM
rkO <- apply(-obsM, 2, frank); rkS <- apply(-stdM, 2, frank); rownames(rkO) <- rownames(rkS) <- rownames(obsM)
pt_rank_o <- frank(-P$st_obs); pt_rank_s <- frank(-P$st_std)
n_st <- a[, .N, by = hv024]
t3 <- data.table(hv024 = as.integer(names(P$st_obs)), state = state_lab[names(P$st_obs)],
  obs = P$st_obs, obs_lo = apply(obsM, 1, qq)[1, ], obs_hi = apply(obsM, 1, qq)[2, ],
  std = P$st_std, std_lo = apply(stdM, 1, qq)[1, ], std_hi = apply(stdM, 1, qq)[2, ],
  change = P$st_std - P$st_obs, ch_lo = apply(chM, 1, qq)[1, ], ch_hi = apply(chM, 1, qq)[2, ],
  rank_obs = pt_rank_o, rank_obs_lo = round(apply(rkO, 1, qq)[1, ]), rank_obs_hi = round(apply(rkO, 1, qq)[2, ]),
  rank_std = pt_rank_s, rank_std_lo = round(apply(rkS, 1, qq)[1, ]), rank_std_hi = round(apply(rkS, 1, qq)[2, ]),
  p_rank_changed = rowMeans(rkO != rkS))
t3 <- merge(t3, n_st, by = "hv024")
t3[, ch_excludes_0 := ch_lo > 0 | ch_hi < 0]
fwrite(t3[order(rank_obs)], file.path(OUT_DIR, "table3_bootstrap.csv"))

pl <- do.call(rbind, lapply(R, `[[`, "pooled")); bs <- sapply(R, `[[`, "bS"); bd <- sapply(R, `[[`, "bD")
po <- data.table(quantity = c("obs", "std", "change", "beta_S", "beta_D"),
  point = c(P$pooled, P$pooled[2] - P$pooled[1], P$bS, P$bD),
  lo = c(qq(pl[, 1])[1], qq(pl[, 2])[1], qq(pl[, 2] - pl[, 1])[1], qq(bs)[1], qq(bd)[1]),
  hi = c(qq(pl[, 1])[2], qq(pl[, 2])[2], qq(pl[, 2] - pl[, 1])[2], qq(bs)[2], qq(bd)[2]))
fwrite(po, file.path(OUT_DIR, "boot_pooled.csv")); print(po)

reg("BOOT_B", length(R))
reg("POOLED_OBS_CI", paste0(f1(po$lo[1]), "\u2013", f1(po$hi[1])))
reg("POOLED_STD_CI", paste0(f1(po$lo[2]), "\u2013", f1(po$hi[2])))
reg("POOLED_CH_CI", ci(po$point[3], po$lo[3], po$hi[3], 1))
reg("BOOT_BS_CI", ci(po$point[4], po$lo[4], po$hi[4]))
for (s in c("Punjab", "Delhi", "Uttarakhand", "Rajasthan", "Uttar Pradesh", "Himachal Pradesh",
            "Mizoram", "West Bengal", "Assam", "Meghalaya")) {
  z <- t3[state == s]; k <- toupper(gsub("[^A-Za-z]", "", s))
  reg(paste0("CH_", k), ci(z$change, z$ch_lo, z$ch_hi, 1))
  rk <- function(p, l, h) sprintf("%d (%d\u2013%d)", as.integer(round(p)), as.integer(l), as.integer(h))
  reg(paste0("RANKO_", k), rk(z$rank_obs, z$rank_obs_lo, z$rank_obs_hi))
  reg(paste0("RANKS_", k), rk(z$rank_std, z$rank_std_lo, z$rank_std_hi))
}
reg("N_CH_SIG", t3[state != "Ladakh" & ch_excludes_0 == TRUE, .N])
reg("N_CH_SIG_NEG", t3[state != "Ladakh" & ch_excludes_0 == TRUE & change < 0, .N])
reg("N_CH_SIG_POS", t3[state != "Ladakh" & ch_excludes_0 == TRUE & change > 0, .N])
reg("RANK_CI_MEDIAN_WIDTH", median(t3$rank_std_hi - t3$rank_std_lo))
reg("N_RANK_CHANGED_90", t3[p_rank_changed >= 0.9, .N])
log_msg("06 done")
