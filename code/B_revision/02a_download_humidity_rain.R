# =============================================================================
# 02a_download_humidity_rain.R — ERA5-Land dewpoint and precipitation
# Reviewer comments 15 and 19 (time-varying confounding by humidity and rain)
#
# Uses the same route that worked for 2-m temperature (ERA_download_fullset.R):
# the hourly dataset "reanalysis-era5-land" (faster queue than daily statistics),
# converted to India local days (UTC+05:00) with the 5-hour "tail" of the
# previous day. Run once; it skips files that already exist.
#   climate/d2m_YYYY_MM.nc  daily mean 2-m dewpoint temperature (K), local days
#   climate/tp_YYYY_MM.nc   daily total precipitation (m), from 00 UTC accumulations
# ERA5-Land precipitation is accumulated from 00 UTC, so the 00 UTC value on day
# t+1 is the total for UTC day t (05:30-05:30 IST); documented in Suppl. Methods S2.
# =============================================================================
source("code/B_revision/00_config.R")
MONTHS <- format(seq(as.Date("2019-05-01"), as.Date("2021-06-01"), by = "month"), "%Y_%m")
need <- c(file.path(CLIM_DIR, paste0("d2m_", MONTHS, ".nc")), file.path(CLIM_DIR, paste0("tp_", MONTHS, ".nc")))
if (all(file.exists(need))) { log_msg("All dewpoint and precipitation files present - nothing to download"); quit(save = "no", status = 0) }
library(terra); library(ecmwfr)
wf_set_key(key = Sys.getenv("CDS_KEY", "PASTE_YOUR_CDS_API_TOKEN"))   # never commit your key
AREA <- c(38, 68, 6, 98)
first_day <- function(ym) as.Date(paste0(sub("_", "-", ym), "-01"))
get <- function(req) {
  if (file.exists(file.path(CLIM_DIR, req$target))) return(invisible())
  wf_request(req, path = CLIM_DIR, transfer = TRUE, time_out = 12 * 3600)
}

## ---- dewpoint: hourly -> local-day mean -------------------------------------
for (ym in MONTHS) {
  out <- file.path(CLIM_DIR, paste0("d2m_", ym, ".nc"))
  if (file.exists(out)) next
  prev <- first_day(ym) - 1
  tryCatch({
    get(list(dataset_short_name = "reanalysis-era5-land", variable = "2m_dewpoint_temperature",
             year = format(prev, "%Y"), month = format(prev, "%m"), day = format(prev, "%d"),
             time = sprintf("%02d:00", 19:23), data_format = "netcdf",
             download_format = "unarchived", area = AREA, target = paste0("d2m_tail_", ym, ".nc")))
    get(list(dataset_short_name = "reanalysis-era5-land", variable = "2m_dewpoint_temperature",
             year = substr(ym, 1, 4), month = substr(ym, 6, 7), day = sprintf("%02d", 1:31),
             time = sprintf("%02d:00", 0:23), data_format = "netcdf",
             download_format = "unarchived", area = AREA, target = paste0("d2m_hourly_", ym, ".nc")))
    r <- c(rast(file.path(CLIM_DIR, paste0("d2m_tail_", ym, ".nc"))),
           rast(file.path(CLIM_DIR, paste0("d2m_hourly_", ym, ".nc"))))
    day <- as.Date(format(time(r) + 5 * 3600, "%Y-%m-%d", tz = "UTC"))
    k <- format(day, "%Y_%m") == ym
    r <- r[[which(k)]]; day <- day[k]; ud <- unique(day)
    dm <- tapp(r, index = match(day, ud), fun = mean); time(dm) <- ud
    writeCDF(dm, out, varname = "d2m", unit = "K", compression = 4, overwrite = TRUE)
    file.remove(file.path(CLIM_DIR, c(paste0("d2m_hourly_", ym, ".nc"), paste0("d2m_tail_", ym, ".nc"))))
    log_msg(ym, "dewpoint done")
  }, error = function(e) log_msg(ym, "dewpoint FAILED:", conditionMessage(e)))
}

## ---- precipitation: 00 UTC accumulation = previous UTC day total ------------
for (ym in MONTHS) {
  out <- file.path(CLIM_DIR, paste0("tp_", ym, ".nc"))
  if (file.exists(out)) next
  nxt <- seq(first_day(ym), by = "month", length.out = 2)[2]
  tryCatch({
    get(list(dataset_short_name = "reanalysis-era5-land", variable = "total_precipitation",
             year = substr(ym, 1, 4), month = substr(ym, 6, 7), day = sprintf("%02d", 2:31),
             time = "00:00", data_format = "netcdf", download_format = "unarchived",
             area = AREA, target = paste0("tp00_", ym, ".nc")))
    get(list(dataset_short_name = "reanalysis-era5-land", variable = "total_precipitation",
             year = format(nxt, "%Y"), month = format(nxt, "%m"), day = "01",
             time = "00:00", data_format = "netcdf", download_format = "unarchived",
             area = AREA, target = paste0("tp00_first_", ym, ".nc")))
    r <- c(rast(file.path(CLIM_DIR, paste0("tp00_", ym, ".nc"))),
           rast(file.path(CLIM_DIR, paste0("tp00_first_", ym, ".nc"))))
    time(r) <- as.Date(time(r)) - 1          # label with the day the total refers to
    writeCDF(r, out, varname = "tp", unit = "m", compression = 4, overwrite = TRUE)
    log_msg(ym, "precipitation done")
  }, error = function(e) log_msg(ym, "precipitation FAILED:", conditionMessage(e)))
}
print(data.frame(month = MONTHS,
                 dewpoint = file.exists(file.path(CLIM_DIR, paste0("d2m_", MONTHS, ".nc"))),
                 precip = file.exists(file.path(CLIM_DIR, paste0("tp_", MONTHS, ".nc")))))
