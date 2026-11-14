# =============================================================================
# 07_lisa.R — district hotspots before/after standardisation, with bootstrap
#             reproducibility of every classification (comments 27, 28)
# Neighbours: 6 nearest district centroids, row-standardised; local Moran's I
# with 999 conditional permutations; P < 0.05; quadrant by mean (spdep).
# Inputs : district_prevalence.csv (05), boot_reps.rds (06), DIST_SHP
# Outputs: lisa_point.csv, lisa_bootstrap_by_district.csv, lisa_summary.csv,
#          figures/FigS6_lisa_maps.png
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages({ library(sf); library(spdep); library(ggplot2) })

dp <- fread(file.path(OUT_DIR, "district_prevalence.csv"))
dp <- dp[n >= SET$min_district_n & hv024 != 37]                       # Ladakh excluded
shp <- st_read(DIST_SHP, quiet = TRUE)
shp$shdist <- as.integer(as.character(shp[[DIST_ID_FIELD]]))
shp <- shp[shp$shdist %in% dp$shdist, ]
dp <- dp[shdist %in% shp$shdist][order(shdist)]
shp <- shp[match(dp$shdist, shp$shdist), ]
reg("LISA_N_DIST", nrow(dp))

sf_use_s2(FALSE)
xy <- st_coordinates(st_point_on_surface(st_geometry(shp)))
lw <- nb2listw(knn2nb(knearneigh(xy, k = SET$lisa_k)), style = "W")

lisa_class <- function(x, nsim = SET$lisa_nsim) {
  lm <- localmoran_perm(x, lw, nsim = nsim)
  p <- lm[, grep("Sim", colnames(lm))[1]]
  q <- as.character(attr(lm, "quadr")$mean)
  cls <- ifelse(p < SET$lisa_alpha, q, "Not significant")
  c("High-High" = "Hotspot", "Low-Low" = "Coldspot", "High-Low" = "High-Low outlier",
    "Low-High" = "Low-High outlier", "Not significant" = "Not significant")[cls]
}
set.seed(SET$seed)
dp[, `:=`(class_obs = lisa_class(obs), class_std = lisa_class(std))]
dp[, changed := class_obs != class_std]
gm <- c(obs = moran.mc(dp$obs, lw, 999)$statistic, std = moran.mc(dp$std, lw, 999)$statistic,
        change = moran.mc(dp$std - dp$obs, lw, 999)$statistic)
reg("MORAN_OBS", f2(gm[1])); reg("MORAN_STD", f2(gm[2])); reg("MORAN_CH", f2(gm[3]))
reg("LISA_CHANGED", dp[changed == TRUE, .N]); reg("LISA_CHANGED_PCT", f1(100 * dp[, mean(changed)]))
reg("HOT_OBS", dp[class_obs == "Hotspot", .N])
reg("HOT_LOST", dp[class_obs == "Hotspot" & class_std != "Hotspot", .N])
reg("HOT_NEW", dp[class_obs != "Hotspot" & class_std == "Hotspot", .N])
reg("DELHI_HOT_OBS", dp[hv024 == 7 & class_obs == "Hotspot", .N])
reg("DELHI_HOT_LOST", dp[hv024 == 7 & class_obs == "Hotspot" & class_std != "Hotspot", .N])
fwrite(dp, file.path(OUT_DIR, "lisa_point.csv"))

## ---- bootstrap reproducibility -------------------------------------------------
bf <- file.path(OUT_DIR, "boot_reps.rds")
if (file.exists(bf)) {
  B <- readRDS(bf)$reps
  key <- as.character(dp$shdist)
  same_o <- same_s <- chg <- matrix(NA, nrow(dp), length(B))
  delhi_all_lost <- n_changed <- numeric(length(B))
  delhi_hot <- dp[, hv024 == 7 & class_obs == "Hotspot"]
  for (b in seq_along(B)) {
    set.seed(SET$seed + b)
    co <- lisa_class(B[[b]]$d_obs[key]); cs <- lisa_class(B[[b]]$d_std[key])
    same_o[, b] <- co == dp$class_obs; same_s[, b] <- cs == dp$class_std; chg[, b] <- co != cs
    n_changed[b] <- sum(co != cs)
    delhi_all_lost[b] <- all(cs[delhi_hot] != "Hotspot")
    if (b %% 50 == 0) log_msg("LISA replicate", b)
  }
  dp[, `:=`(pct_same_obs = 100 * rowMeans(same_o), pct_same_std = 100 * rowMeans(same_s),
            pct_reps_changed = 100 * rowMeans(chg))]
  fwrite(dp, file.path(OUT_DIR, "lisa_bootstrap_by_district.csv"))
  ls <- data.table(
    replicates = length(B),
    median_n_changed = median(n_changed), n_changed_lo = quantile(n_changed, .025),
    n_changed_hi = quantile(n_changed, .975),
    reclassified_in_50pct = dp[changed == TRUE & pct_reps_changed >= 50, .N],
    reclassified_in_80pct = dp[changed == TRUE & pct_reps_changed >= 80, .N],
    hotspot_obs_reproduced_80 = dp[class_obs == "Hotspot" & pct_same_obs >= 80, .N],
    pct_reps_all_delhi_lost = 100 * mean(delhi_all_lost))
  fwrite(ls, file.path(OUT_DIR, "lisa_summary.csv")); print(ls)
  reg("LISA_BOOT_CH_CI", sprintf("%d (95%% interval %d\u2013%d)", as.integer(ls$median_n_changed),
                                 as.integer(ls$n_changed_lo), as.integer(ls$n_changed_hi)))
  reg("LISA_REPRO50", ls$reclassified_in_50pct); reg("LISA_REPRO80", ls$reclassified_in_80pct)
  reg("DELHI_LOST_PCT", f1(ls$pct_reps_all_delhi_lost))
}

## ---- basic maps ------------------------------------------------------------------
mk <- function(lbl, cls) { m <- shp[, "shdist"]; m$panel <- lbl; m$cls <- cls; m }
mp <- rbind(mk("observed", dp$class_obs), mk("standardised", dp$class_std))
ggsave(file.path(FIG_DIR, "lisa_classes.png"), ggplot(mp) + geom_sf(aes(fill = cls)) + facet_wrap(~panel))
log_msg("07 done")
