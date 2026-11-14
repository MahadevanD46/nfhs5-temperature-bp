###############################################################################
# ERA5-Land: download hourly 2 m temperature and convert to daily means
# (India local days, UTC+05:00). Run the whole script top to bottom.
###############################################################################
library(terra)
library(ecmwfr)

## ===================== SETTINGS: edit these two lines =======================
PROJECT_DIR <- "PASTE_YOUR_PROJECT_FOLDER"            # <- paste the result of getwd()
CDS_KEY     <- "PASTE_YOUR_CDS_API_TOKEN"     # <- your CDS token
## ============================================================================

MONTHS <- c("2019_07", "2019_08", "2019_09", "2019_10", "2019_11",
            "2019_12", "2020_11", "2020_12", "2021_05")

CLIM_DIR <- file.path(PROJECT_DIR, "climate")
if (!dir.exists(CLIM_DIR))
  stop("Folder not found: ", CLIM_DIR, "  -> check PROJECT_DIR")
if (!file.exists(file.path(CLIM_DIR, "t2m_monthly_1991_2018.nc")))
  stop("This is not your climate folder (t2m_monthly_1991_2018.nc not found): ", CLIM_DIR)

wf_set_key(key = CDS_KEY)

# Delete half-finished temporary files from earlier failed attempts
unlink(list.files(CLIM_DIR, pattern = "^ecmwfr_", full.names = TRUE))

## ---------------------------- helper functions ------------------------------
f_daily  <- function(ym) file.path(CLIM_DIR, paste0("t2m_", ym, ".nc"))
f_hourly <- function(ym) file.path(CLIM_DIR, paste0("t2m_hourly_", ym, ".nc"))
f_tail   <- function(ym) file.path(CLIM_DIR, paste0("t2m_tail_", ym, ".nc"))
first_day <- function(ym) as.Date(paste0(sub("_", "-", ym), "-01"))
next_ym   <- function(ym) format(seq(first_day(ym), by = "month", length.out = 2)[2], "%Y_%m")

# Full month of hourly data
req_hourly <- function(ym) list(
  dataset_short_name = "reanalysis-era5-land",
  variable = "2m_temperature",
  year  = substr(ym, 1, 4),
  month = substr(ym, 6, 7),
  day   = sprintf("%02d", 1:31),
  time  = sprintf("%02d:00", 0:23),
  data_format = "netcdf", download_format = "unarchived",
  area = c(38, 68, 6, 98),
  target = basename(f_hourly(ym)))

# Small request: 19:00-23:00 UTC of the previous month's last day
# (= 00:00-04:00 India time on day 1 of this month)
req_tail <- function(ym) {
  prev <- first_day(ym) - 1
  list(
    dataset_short_name = "reanalysis-era5-land",
    variable = "2m_temperature",
    year  = format(prev, "%Y"),
    month = format(prev, "%m"),
    day   = format(prev, "%d"),
    time  = sprintf("%02d:00", 19:23),
    data_format = "netcdf", download_format = "unarchived",
    area = c(38, 68, 6, 98),
    target = basename(f_tail(ym)))
}

download <- function(req) {
  wf_request(req, path = CLIM_DIR, transfer = TRUE, time_out = 12 * 3600)
  if (!file.exists(file.path(CLIM_DIR, req$target)))
    stop("Download finished but file not found: ", req$target)
}

# Save the last 5 hours of a month as the 'tail' for the next month
save_tail_for_next <- function(r, ym) {
  ft <- f_tail(next_ym(ym))
  if (!file.exists(ft))
    writeCDF(r[[(nlyr(r) - 4):nlyr(r)]], ft, varname = "t2m", unit = "K", overwrite = TRUE)
}

# Check that an existing tail file is usable
tail_ok <- function(ym) {
  ft <- f_tail(ym)
  if (!file.exists(ft)) return(FALSE)
  tt <- time(rast(ft))
  good <- length(tt) == 5 && !anyNA(tt)
  if (!good) file.remove(ft)
  good
}

# Hourly (UTC) -> daily mean over India local days
hourly_to_daily <- function(ym) {
  r <- rast(f_hourly(ym))
  save_tail_for_next(r, ym)
  r <- c(rast(f_tail(ym)), r)
  day  <- as.Date(format(time(r) + 5 * 3600, "%Y-%m-%d", tz = "UTC"))
  keep <- format(day, "%Y_%m") == ym
  r <- r[[which(keep)]]; day <- day[keep]
  ud  <- unique(day)
  hrs <- tabulate(match(day, ud))
  dm  <- tapp(r, index = match(day, ud), fun = mean)
  time(dm) <- ud
  writeCDF(dm, f_daily(ym), varname = "t2m", unit = "K", compression = 4, overwrite = TRUE)
  cat(ym, ": ", length(ud), " days; hours per day = ",
      paste(unique(hrs), collapse = "/"), "\n", sep = "")
}

## ------------------------------- run ----------------------------------------
# Re-use the June 2019 hourly file you already have as the tail for July 2019
if (file.exists(f_hourly("2019_06")) && !tail_ok("2019_07"))
  save_tail_for_next(rast(f_hourly("2019_06")), "2019_06")

for (ym in MONTHS) {
  if (file.exists(f_daily(ym))) { cat(ym, "already done\n"); next }
  tryCatch({
    if (!tail_ok(ym))               download(req_tail(ym))
    if (!file.exists(f_hourly(ym))) download(req_hourly(ym))
    hourly_to_daily(ym)
    file.remove(f_hourly(ym), f_tail(ym))          # free disk space
  }, error = function(e) cat(ym, "FAILED:", conditionMessage(e), "\n"))
}

## ------------------------------ summary -------------------------------------
print(data.frame(month = MONTHS, done = file.exists(f_daily(MONTHS))))
