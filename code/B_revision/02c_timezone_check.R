# =============================================================================
# 02c_timezone_check.R — hourly -> daily computation and time zone (comment 11)
#
# ERA5-Land 2-m temperature is an instantaneous field on whole UTC hours.
# A local day in IST (UTC+05:30) runs 18:30 UTC to 18:30 UTC; the hourly
# instants inside it are 19:00, 20:00, ..., 18:00 UTC — exactly the 24
# instants inside a UTC+05:00 day (19:00 to 19:00 UTC). The daily mean from
# hourly data is therefore IDENTICAL under UTC+05:00 and UTC+05:30. This script
# demonstrates that numerically, shows how much a one-hour shift (UTC+06:00) or
# a UTC day would change the mean, and compares the CDS daily-statistics
# product with the mean computed from hourly data.
# Needs: climate/t2m_hourly_2020_01.nc  (+ tail file 2019-12-31 19:00-23:00 UTC)
#        climate/t2m_2020_01.nc         (daily-statistics product used in the paper)
# Output: output_revision/timezone_check.csv -> Supplementary Table S6
# =============================================================================
source("code/B_revision/00_config.R")
library(terra)
fh <- file.path(CLIM_DIR, "t2m_hourly_2020_01.nc"); ft <- file.path(CLIM_DIR, "t2m_tail_2020_01.nc")
if (!file.exists(fh)) {
  library(ecmwfr); wf_set_key(key = Sys.getenv("CDS_KEY", "PASTE_YOUR_CDS_API_TOKEN"))
  base <- list(dataset_short_name = "reanalysis-era5-land", variable = "2m_temperature",
               data_format = "netcdf", download_format = "unarchived", area = c(38, 68, 6, 98))
  wf_request(c(base, list(year = "2020", month = "01", day = sprintf("%02d", 1:31),
               time = sprintf("%02d:00", 0:23), target = basename(fh))),
             path = CLIM_DIR, transfer = TRUE, time_out = 12 * 3600)
  wf_request(c(base, list(year = "2019", month = "12", day = "31",
               time = sprintf("%02d:00", 18:23), target = basename(ft))),
             path = CLIM_DIR, transfer = TRUE, time_out = 12 * 3600)
}
r <- c(rast(ft), rast(fh)); tt <- time(r)
land <- which(!is.na(values(r[[1]], mat = FALSE)))
cells <- sort(sample(land, min(5000, length(land))))
H <- as.matrix(terra::extract(r, cells)) - 273.15

daily <- function(offset_h) {
  day <- as.Date(format(tt + offset_h * 3600, "%Y-%m-%d", tz = "UTC"))
  full <- names(which(table(day) == 24))
  full <- full[substr(full, 1, 7) == "2020-01"]
  sapply(full, function(d) rowMeans(H[, day == as.Date(d), drop = FALSE]))
}
d500 <- daily(5); d530 <- daily(5.5); d600 <- daily(6); d000 <- daily(0)
prod <- rast(file.path(CLIM_DIR, "t2m_2020_01.nc"))
P <- as.matrix(terra::extract(prod, cells)) - 273.15
colnames(P) <- as.character(as.Date(time(prod)))
days <- Reduce(intersect, list(colnames(d500), colnames(d530), colnames(d600), colnames(d000), colnames(P)))
cmp <- function(lbl, X) data.table(comparison = lbl,
  mean_abs_diff = mean(abs(X[, days] - d500[, days])), max_abs_diff = max(abs(X[, days] - d500[, days])),
  correlation = cor(c(X[, days]), c(d500[, days])), cell_days = length(cells) * length(days))
tz <- rbind(cmp("Hourly, IST day (UTC+05:30) vs UTC+05:00 day", d530),
            cmp("Hourly, UTC+06:00 day vs UTC+05:00 day", d600),
            cmp("Hourly, UTC day vs UTC+05:00 day", d000),
            cmp("CDS daily-statistics product (UTC+05:00) vs hourly UTC+05:00", P))
fwrite(tz, file.path(OUT_DIR, "timezone_check.csv")); print(tz)
reg("TZ_MAD_IST", formatC(tz$mean_abs_diff[1], format = "f", digits = 3))
reg("TZ_MAD_1H", formatC(tz$mean_abs_diff[2], format = "f", digits = 2))
reg("TZ_MAD_PROD", formatC(tz$mean_abs_diff[4], format = "f", digits = 3))
reg("TZ_COR_PROD", formatC(tz$correlation[4], format = "f", digits = 4))
