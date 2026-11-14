# Synthetic NFHS-5-like project for testing the pipeline. NOT real data: all values are simulated.
# Creates synthetic_project/ with IAPR7EFL.DTA, IAGE7AFL.shp, climate/*.nc and boundaries/ using the
# same file and variable names as the real inputs. Planted effects: SBP -0.23, DBP -0.12 mmHg per degC.
suppressPackageStartupMessages({ library(data.table); library(terra); library(sf); library(haven) })
set.seed(42)
dir.create("synthetic_project", showWarnings = FALSE); setwd("synthetic_project")
dir.create("climate", showWarnings = FALSE); dir.create("boundaries", showWarnings = FALSE)

## ---- grid and land mask (0.1 deg, India box) --------------------------------
tmpl <- rast(xmin = 68, xmax = 98, ymin = 6, ymax = 38, resolution = 0.1, crs = "EPSG:4326")
xy <- xyFromCell(tmpl, 1:ncell(tmpl))
land <- !( (xy[,1] > 88 & xy[,2] < 20) | (xy[,1] < 72 & xy[,2] < 22) | (xy[,2] < 8) |
           (xy[,1] > 80 & xy[,1] < 86 & xy[,2] < 14) )
lat <- xy[, 2]; lon <- xy[, 1]
base_annual <- 27 - 0.45 * pmax(lat - 15, 0) - 8 * (lat > 32 & lon < 80)     # cold north-west mountains
amp <- 2 + 0.35 * pmax(lat - 10, 0)                                           # seasonal amplitude
mk <- function(dates, noise_sd = 2.5) {
  doy <- as.integer(format(dates, "%j"))
  r <- rast(tmpl, nlyrs = length(dates))
  seas <- -cos(2 * pi * (doy - 15) / 365.25)                                  # coldest mid-Jan
  ar <- as.numeric(arima.sim(list(ar = .8), length(dates))) * noise_sd * .6
  M <- outer(base_annual, rep(1, length(dates))) + outer(amp, seas) * 1 +
       matrix(rep(ar, each = ncell(tmpl)), ncell(tmpl)) + rnorm(ncell(tmpl) * length(dates), 0, .7)
  M[!land, ] <- NA
  values(r) <- M + 273.15; time(r) <- dates; r
}
months <- seq(as.Date("2019-05-01"), as.Date("2021-06-01"), by = "month")
CLIM_DONE <- file.exists("climate/t2m_tail_2020_01.nc")
if (!CLIM_DONE) for (i in seq_along(months)) {
  m <- months[i]; dd <- seq(m, seq(m, by = "month", length.out = 2)[2] - 1, by = "day")
  r <- mk(dd)
  writeCDF(r, sprintf("climate/t2m_%s.nc", format(m, "%Y_%m")), varname = "t2m", unit = "K", overwrite = TRUE)
  dp <- r - runif(1, 3, 12); time(dp) <- dd; writeCDF(dp, sprintf("climate/d2m_%s.nc", format(m, "%Y_%m")), varname = "d2m", unit = "K", overwrite = TRUE)
  tp <- rast(r); vv <- matrix(rexp(ncell(r) * nlyr(r), 400), ncell(r)); vv[!land, ] <- NA; values(tp) <- vv; time(tp) <- dd
  writeCDF(tp, sprintf("climate/tp_%s.nc", format(m, "%Y_%m")), varname = "tp", unit = "m", overwrite = TRUE)
}
# monthly normals 1991-2018 (336 layers)
if (!CLIM_DONE) {
mm <- seq(as.Date("1991-01-15"), as.Date("2018-12-15"), by = "month")
rn <- rast(tmpl, nlyrs = length(mm))
seas <- -cos(2 * pi * (as.integer(format(mm, "%j")) - 15) / 365.25)
M <- outer(base_annual, rep(1, length(mm))) + outer(amp, seas) + rnorm(ncell(tmpl) * length(mm), 0, .3)
M[!land, ] <- NA; values(rn) <- M + 273.15; time(rn) <- as.Date(format(mm, "%Y-%m-01"))
writeCDF(rn, "climate/t2m_monthly_1991_2018.nc", varname = "t2m", unit = "K", overwrite = TRUE)
# hourly January 2020 + tail
hh <- seq(as.POSIXct("2020-01-01 00:00", tz = "UTC"), as.POSIXct("2020-01-31 23:00", tz = "UTC"), by = "hour")
mkh <- function(tt) { r <- rast(tmpl, nlyrs = length(tt)); hr <- as.integer(format(tt, "%H"))
  b <- base_annual - amp + 273.15; b[!land] <- NA
  M <- matrix(b, ncell(tmpl), length(tt)); M <- M + rep(4 * sin(2 * pi * (hr - 3) / 24), each = ncell(tmpl))
  values(r) <- M; rm(M); gc(); time(r) <- tt; r }
tmpf <- character()
for (dd in split(hh, as.Date(hh))) { f <- tempfile(fileext = ".tif"); writeRaster(mkh(dd), f); tmpf <- c(tmpf, f) }
rh <- rast(tmpf); time(rh) <- hh
writeCDF(rh, "climate/t2m_hourly_2020_01.nc", varname = "t2m", unit = "K", overwrite = TRUE)
writeCDF(mkh(seq(as.POSIXct("2019-12-31 18:00", tz = "UTC"), by = "hour", length.out = 6)),
         "climate/t2m_tail_2020_01.nc", varname = "t2m", unit = "K", overwrite = TRUE)
}

## ---- districts (2-degree squares) and states ------------------------------------
sq <- st_make_grid(st_as_sfc(st_bbox(c(xmin = 68, xmax = 98, ymin = 8, ymax = 38), crs = 4326)), cellsize = 2)
sq <- st_sf(geometry = sq); cen <- st_coordinates(st_centroid(sq))
landfrac <- sapply(seq_len(nrow(sq)), function(i) mean(land[cellFromXY(tmpl, cbind(runif(50, cen[i,1]-1, cen[i,1]+1), runif(50, cen[i,2]-1, cen[i,2]+1)))]))
sq <- sq[landfrac > .5, ]; cen <- cen[landfrac > .5, ]
sq$REGCODE <- seq_len(nrow(sq))
st_codes <- c(1:25, 27:30, 32:37)                       # 35 states (31 Lakshadweep handled apart)
ord <- order(-cen[, 2], cen[, 1])
sq$hv024 <- NA; sq$hv024[ord] <- st_codes[pmin(length(st_codes), ceiling(seq_along(ord) / (length(ord) / length(st_codes))))]
st_write(sq["REGCODE"], "boundaries/sdr_subnational_boundaries2.shp", quiet = TRUE, delete_dsn = TRUE)

## ---- clusters --------------------------------------------------------------------
ncl_d <- 8
cl <- rbindlist(lapply(seq_len(nrow(sq)), function(i) data.table(shdist = sq$REGCODE[i], hv024 = sq$hv024[i],
  LONGNUM = runif(ncl_d, cen[i,1] - .95, cen[i,1] + .95), LATNUM = runif(ncl_d, cen[i,2] - .95, cen[i,2] + .95))))
cl[, on_land := land[cellFromXY(tmpl, cbind(LONGNUM, LATNUM))]]
cl <- cl[on_land == TRUE]
# move a few clusters just off the coast (to test 3x3 / 5x5 infill)
coast <- sample(which(cl$LONGNUM > 86 & cl$LATNUM < 21 & cl$LATNUM > 19), 4)
cl[coast, LONGNUM := 88.05]
cl[coast[1:2], LONGNUM := 88.25]
# Lakshadweep: sea, unfillable
cl <- rbind(cl, data.table(shdist = 999, hv024 = 31, LONGNUM = 72.6, LATNUM = 10.5, on_land = FALSE), fill = TRUE)
cl[, hv001 := .I + 10000]
cl[, URBAN_RURA := sample(c("U", "R"), .N, TRUE, c(.3, .7))]
cl[, ALT_DEM := fifelse(LATNUM > 32, runif(.N, 1500, 4000), runif(.N, 0, 800))]
nogps <- sample(cl$hv001, 3)
ge <- st_as_sf(cl[, .(DHSCLUST = hv001, URBAN_RURA, ALT_DEM, LATNUM = fifelse(hv001 %in% nogps, 0, LATNUM),
                      LONGNUM = fifelse(hv001 %in% nogps, 0, LONGNUM))], coords = c("LONGNUM", "LATNUM"),
               crs = 4326, remove = FALSE)
st_write(ge, "IAGE7AFL.shp", quiet = TRUE, delete_dsn = TRUE)

## ---- dates: phase by state -----------------------------------------------------------
ph2 <- c(1:9, 12, 20, 21:23, 33, 34, 35, 37)
cl[, shphase := fifelse(hv024 %in% ph2, 2L, 1L)]
pick_date <- function(ph) if (ph == 1) as.Date("2019-06-15") + sample(0:220, 1) else
  sample(c(seq(as.Date("2020-01-02"), as.Date("2020-03-20"), "day"), seq(as.Date("2020-11-05"), as.Date("2021-04-28"), "day")), 1)
cl[, cdate := as.Date(sapply(shphase, pick_date), origin = "1970-01-01")]

## ---- persons ----------------------------------------------------------------------
pp <- cl[, .(hvidx = 1:sample(25:45, 1)), by = .(hv001, hv024, shdist, shphase, cdate, URBAN_RURA, ALT_DEM, LONGNUM, LATNUM, on_land)]
pp[, hv002 := (hvidx - 1) %/% 5 + 1]
N <- nrow(pp)
pp[, `:=`(hv104 = sample(1:2, N, TRUE), hv105 = sample(c(0:14, 15:90), N, TRUE, c(rep(.5, 15), rep(1.3, 76))),
          hv005 = round(runif(N, 2e5, 2e6)), hv021 = hv001, hv025 = fifelse(URBAN_RURA == "U", 1L, 2L),
          hv270 = sample(1:5, N, TRUE), hv040 = round(ALT_DEM), hv018 = sample(1:6, N, TRUE),
          hv030 = 1L, hv027 = 1L, hv206 = 1L, hv215 = 1L, hv103 = sample(c(1L, 0L), N, TRUE, c(.95, .05)))]
pp[, hv022 := shdist * 10 + hv025]
pp[, bio := cdate + sample(c(0, 0, 1, 1, 2, 3, 5, 9), N, TRUE)]
pp[, int := cdate]
pp[sample(N, 25), bio := as.Date("2020-04-10")]            # suspension "date errors"
pp[, `:=`(hv807d = as.integer(format(bio, "%d")), hv807m = as.integer(format(bio, "%m")), hv807y = as.integer(format(bio, "%Y")),
          hv016 = as.integer(format(int, "%d")), hv006 = as.integer(format(int, "%m")), hv007 = as.integer(format(int, "%Y")))]
pp[sample(N, 150), `:=`(hv807d = NA, hv807m = NA, hv807y = NA)]
# temperature on bio day for the outcome model
cell <- cellFromXY(tmpl, cbind(pp$LONGNUM, pp$LATNUM))
bd <- fcoalesce(as.Date(sprintf("%04d-%02d-%02d", pp$hv807y, pp$hv807m, pp$hv807d), optional = TRUE), pp$int)
tv <- rep(NA_real_, N)
for (mi in seq_along(months)) { m <- months[mi]
  r <- rast(sprintf("climate/t2m_%s.nc", format(m, "%Y_%m"))); dts <- as.Date(time(r))
  k <- which(bd %in% dts & !is.na(cell))
  if (!length(k)) next
  v <- terra::extract(r, cell[k])
  tv[k] <- as.matrix(v)[cbind(seq_along(k), match(bd[k], dts))] - 273.15
}
tv[is.na(tv)] <- 25
age <- pmax(pp$hv105, 15)
sbp_true <- 105 + .45 * (age - 15) + 4 * (pp$hv104 == 1) - 0.23 * tv + rnorm(N, 0, 15)
dbp_true <- 70 + .15 * (age - 15) + 2 * (pp$hv104 == 1) - 0.12 * tv + .3 * (sbp_true - 120) + rnorm(N, 0, 7)
rd <- function(x) round(x + rnorm(length(x), 0, 3))
pp[, `:=`(shb18s = rd(sbp_true + 3), shb18d = rd(dbp_true + 1), shb25s = rd(sbp_true), shb25d = rd(dbp_true),
          shb29s = rd(sbp_true), shb29d = rd(dbp_true))]
child <- pp$hv105 < 15; nobp <- sample(N, round(.04 * N))
for (v in c("shb18s", "shb18d", "shb25s", "shb25d", "shb29s", "shb29d")) set(pp, which(child | seq_len(N) %in% nobp), v, NA)
pp[sample(N, 30), shb25s := 300]                           # out of range
pp[, `:=`(shb14a = sample(0:1, N, TRUE, c(.8, .2)), shb14b = sample(0:1, N, TRUE, c(.8, .2)),
          shb14c = sample(0:1, N, TRUE, c(.9, .1)), shb14d = 0L, sh25 = sample(0:1, N, TRUE, c(.75, .25)),
          sh26 = sample(0:1, N, TRUE, c(.9, .1)), sh50g = sample(0:1, N, TRUE), sh50q = sample(0:1, N, TRUE, c(.8, .2)),
          shb17 = sample(c(700:1959), N, TRUE), shb24 = 0L, shb28 = 0L, shb20 = sample(0:1, N, TRUE, c(.85, .15)))]
pp[, shb21 := fifelse(shb20 == 1, sample(0:1, N, TRUE), NA_integer_)]
pp[, shb53 := sample(0:12, N, TRUE)][, shb54 := 0L][, shb57 := sample(0:1, N, TRUE, c(.95, .05))][, shb73 := 1000L]
pp[, shb74 := round(exp(log(110) + 0.003 * tv + .02 * (age - 40) / 10 + rnorm(N, 0, .25)))]
pp[sample(N, 300), shb74 := 996]
pp[, `:=`(ha54 = fifelse(hv104 == 2 & hv105 < 45, sample(0:1, N, TRUE, c(.96, .04)), 0L),
          ha3 = fifelse(hv104 == 2 & hv105 <= 49, round(rnorm(N, 1520, 60)), NA_real_),
          hb3 = fifelse(hv104 == 1 & hv105 <= 54, round(rnorm(N, 1650, 70)), NA_real_))]
pp[, `:=`(ha40 = fifelse(!is.na(ha3), round(rnorm(N, 2200, 300)), NA_real_), hb40 = fifelse(!is.na(hb3), round(rnorm(N, 2200, 300)), NA_real_))]
pp[sample(N, 20), sh25 := NA]                               # missing covariates
keep <- c("hv001","hv002","hvidx","hv005","hv103","hv006","hv007","hv016","hv021","hv022","hv024",
          "hv025","hv027","hv040","hv104","hv105","hv206","hv215","hv270","hv018","hv030",
          "hv807d","hv807m","hv807y","shphase","shdist","sh25","sh26","sh50g","sh50q",
          "shb14a","shb14b","shb14c","shb14d","shb17","shb18s","shb18d","shb20","shb21",
          "shb24","shb25s","shb25d","shb28","shb29s","shb29d",
          "shb53","shb54","shb57","shb73","shb74","ha3","ha40","ha54","hb3","hb40")
write_dta(as.data.frame(pp[, ..keep]), "IAPR7EFL.DTA")
cat("Persons:", N, "| clusters:", nrow(cl), "| districts:", nrow(sq), "\n")
