# =============================================================================
# 02d_imd_validation.R — OPTIONAL independent validation against IMD gridded
# station-based daily temperature (comment 11: "validate against station data
# or another reanalysis if possible")
#
# Data: India Meteorological Department 1.0 x 1.0 degree gridded daily maximum
# and minimum temperature (Srivastava et al., Atmos Sci Lett 2009), free from
# https://www.imdpune.gov.in (Gridded data -> Temperature). Save the yearly
# binary files as imd/Maxtemp_MaxT_YYYY.GRD and imd/Mintemp_MinT_YYYY.GRD.
# Grid (check the IMD read-me that comes with the files): 31 x 31 cells,
# longitude 67.5-97.5 E, latitude 7.5-37.5 N, float32, missing = 99.9.
# Daily mean is approximated as (Tmax + Tmin) / 2.
# Preferred input: imd/IMD_Tmax_2019_2021.nc and imd/IMD_Tmin_2019_2021.nc written by
# imdlib (Python): imd.get_data('tmax', 2019, 2021).get_xarray().to_netcdf(...).
# Output: output_revision/imd_validation.csv  -> Supplementary Table S6
# =============================================================================
source("code/B_revision/00_config.R")
library(terra)
NX <- 31; NY <- 31; LON0 <- 67.5; LAT0 <- 7.5
read_imd <- function(f, year) {
  nd <- if (year %% 4 == 0) 366 else 365
  v <- readBin(f, "numeric", n = NX * NY * nd, size = 4, endian = "little")
  v[v > 99] <- NA
  arr <- array(v, c(NX, NY, nd))                       # lon fastest, then lat, then day
  r <- rast(nrows = NY, ncols = NX, nlyrs = nd, xmin = LON0 - .5, xmax = LON0 + NX - .5,
            ymin = LAT0 - .5, ymax = LAT0 + NY - .5, crs = "EPSG:4326")
  values(r) <- matrix(arr[, NY:1, ], ncol = nd)     # terra order: north row first, west to east
  time(r) <- seq(as.Date(paste0(year, "-01-01")), by = "day", length.out = nd)
  r
}
if (!dir.exists(IMD_DIR)) { log_msg("IMD files not found - optional validation skipped"); quit(save = "no", status = 0) }

## ---- read IMD: NetCDF from imdlib (preferred) or raw yearly .GRD files ----------
nc_max <- file.path(IMD_DIR, "IMD_Tmax_2019_2021.nc"); nc_min <- file.path(IMD_DIR, "IMD_Tmin_2019_2021.nc")
if (file.exists(nc_max) && file.exists(nc_min)) {
  rd_nc <- function(f) {
    r <- rast(f)
    if (length(varnames(r)) > 1) r <- r[[grep("^t", names(r))]]
    tt <- as.Date(time(r))
    if (anyNA(tt)) stop("No time axis read from ", f, " - send this message to Claude")
    r <- classify(r, cbind(90, Inf, NA))                 # IMD missing value 99.9 -> NA
    time(r) <- tt; r
  }
  tx <- rd_nc(nc_max); tn <- rd_nc(nc_min)
  common <- intersect(as.character(time(tx)), as.character(time(tn)))
  tx <- tx[[match(common, as.character(time(tx)))]]; tn <- tn[[match(common, as.character(time(tn)))]]
  tm <- (tx + tn) / 2; time(tm) <- as.Date(common)
  imd_list <- list(tm)
} else {
  imd_list <- lapply(2019:2021, function(y) {
    tx <- read_imd(file.path(IMD_DIR, sprintf("Maxtemp_MaxT_%d.GRD", y)), y)
    tn <- read_imd(file.path(IMD_DIR, sprintf("Mintemp_MinT_%d.GRD", y)), y)
    tm <- (tx + tn) / 2; time(tm) <- time(tx); tm })
}
# sanity check: Delhi (28.6N, 77.2E) January mean should be roughly 13-16 degC
chk <- imd_list[[1]]; jan <- which(format(as.Date(time(chk)), "%m") == "01")
delhi <- mean(as.numeric(terra::extract(chk[[jan]], cbind(77.2, 28.6))[1, ]), na.rm = TRUE)
log_msg("IMD check: layers =", sum(sapply(imd_list, nlyr)), "| Delhi January mean (Tmax+Tmin)/2 =", round(delhi, 1), "degC")

a <- load_a(cols = c("hv001", "date", "t0", "lat", "lon", "sbp", "dbp", "wt", "age_grp",
                    "hv104", "hv270", "hv025", "shb14a", "shb14b", "shb14c", "sh25", "sh26",
                    "shdist", "ym", "dow"))
cl <- unique(a[, .(hv001, lat, lon)])
imd <- rbindlist(lapply(imd_list, function(tm) {
  v <- as.matrix(terra::extract(tm, cbind(cl$lon, cl$lat)))
  data.table(hv001 = rep(cl$hv001, ncol(v)), date = rep(as.Date(time(tm)), each = nrow(v)), t_imd = as.vector(v))
}))
a <- merge(a, imd, by = c("hv001", "date"))
a <- a[!is.na(t_imd)]
reg("IMD_DELHI_JAN", round(delhi, 1))
m_imd <- fit("sbp", "t_imd", a); m_era <- fit("sbp", "t0", a)
out <- data.table(
  n = nrow(a), correlation = a[, cor(t0, t_imd, use = "complete.obs")],
  mean_diff = a[, mean(t0 - t_imd, na.rm = TRUE)], sd_diff = a[, sd(t0 - t_imd, na.rm = TRUE)],
  beta_imd = coef(m_imd)["t_imd"], se_imd = se(m_imd)["t_imd"],
  beta_era5_same_sample = coef(m_era)["t0"], se_era5 = se(m_era)["t0"])
fwrite(out, file.path(OUT_DIR, "imd_validation.csv")); print(out)
reg("IMD_COR", formatC(out$correlation, format = "f", digits = 2))
reg("IMD_BETA", ci(out$beta_imd, out$beta_imd - 1.96 * out$se_imd, out$beta_imd + 1.96 * out$se_imd))
