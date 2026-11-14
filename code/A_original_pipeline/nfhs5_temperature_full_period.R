###############################################################################
# Measurement-day ambient temperature, blood pressure and glucose in NFHS-5
# FULL FIELDWORK PERIOD (June 2019 - May 2021)
#   Inputs : IAPR7EFL.DTA (NFHS-5 PR file), IAGE7AFL.shp (GPS),
#            climate/t2m_YYYY_MM.nc (daily means, India local days)
#            climate/t2m_monthly_1991_2018.nc (monthly means for normals)
#   Outputs: output_full/ (analytic data, tables, figures, model objects)
# Tip (RStudio): set Environment pane to "Manual Refresh Only".
###############################################################################

## =========================== 0. SETUP =======================================
pkgs <- c("haven", "data.table", "lubridate", "sf", "terra", "fixest", "ggplot2",
          "ncdf4", "Hmisc", "tidyselect")
miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss)
library(haven); library(data.table); library(lubridate)
library(sf); library(terra); library(fixest); library(ggplot2)

PR_FILE  <- "IAPR7EFL.DTA"
GE_FILE  <- "IAGE7AFL.shp"
CLIM_DIR <- "climate"
OUT_DIR  <- "output_full"
dir.create(OUT_DIR, showWarnings = FALSE)

FILL_NORMALS <- TRUE    # coast-fill the normals too (fixes missing t_annual for coastal clusters)

state_lab <- c("1"="Jammu & Kashmir","2"="Himachal Pradesh","3"="Punjab","4"="Chandigarh",
               "5"="Uttarakhand","6"="Haryana","7"="Delhi","8"="Rajasthan","9"="Uttar Pradesh","10"="Bihar",
               "11"="Sikkim","12"="Arunachal Pradesh","13"="Nagaland","14"="Manipur","15"="Mizoram",
               "16"="Tripura","17"="Meghalaya","18"="Assam","19"="West Bengal","20"="Jharkhand","21"="Odisha",
               "22"="Chhattisgarh","23"="Madhya Pradesh","24"="Gujarat","25"="Dadra & Nagar Haveli and Daman & Diu",
               "27"="Maharashtra","28"="Andhra Pradesh","29"="Karnataka","30"="Goa","31"="Lakshadweep",
               "32"="Kerala","33"="Tamil Nadu","34"="Puducherry","35"="Andaman & Nicobar Islands",
               "36"="Telangana","37"="Ladakh")

## ===================== 1. NFHS-5 DATA PREPARATION ==========================
keep <- c("hv001","hv002","hvidx","hv005","hv006","hv007","hv016","hv021","hv022","hv024",
          "hv025","hv027","hv040","hv104","hv105","hv206","hv215","hv270","hv018","hv030",
          "hv807d","hv807m","hv807y","shphase","shdist","sh25","sh26","sh50g","sh50q",
          "shb14a","shb14b","shb14c","shb14d","shb17","shb18s","shb18d","shb20","shb21",
          "shb24","shb25s","shb25d","shb28","shb29s","shb29d",
          "shb53","shb54","shb57","shb73","shb74","ha3","ha40","ha54","hb3","hb40")

pr <- read_dta(PR_FILE, col_select = tidyselect::all_of(keep)) |> zap_labels() |> as.data.table()

ok <- function(x, lo, hi) fifelse(!is.na(x) & x >= lo & x <= hi, as.numeric(x), NA_real_)
hhmm_to_hr <- function(x) { x <- ok(x, 0, 2359); h <- floor(x / 100) + (x %% 100) / 60
fifelse(h >= 6 & h <= 21, h, NA_real_) }

pr[, `:=`(
  s2 = ok(shb25s, 70, 270), s3 = ok(shb29s, 70, 270),
  d2 = ok(shb25d, 40, 150), d3 = ok(shb29d, 40, 150),
  glu     = ok(shb74, 40, 500),
  bmi     = ok(fifelse(hv104 == 2, ha40, hb40), 1200, 6000) / 100,
  height  = ok(fifelse(hv104 == 2, ha3,  hb3),  1000, 2200) / 10,
  bp_med  = as.integer(fcoalesce(as.numeric(shb21), 0) == 1),
  glu_med = as.integer(fcoalesce(as.numeric(shb57), 0) == 1),
  preg    = as.integer(hv104 == 2 & fcoalesce(as.numeric(ha54), 0) == 1),
  bp_hour = hhmm_to_hr(shb17),
  fan     = as.integer(sh50g == 1),
  cooler  = as.integer(sh50q == 1),
  wt      = hv005 / 1e6
)]
pr[, sbp := rowMeans(cbind(s2, s3), na.rm = TRUE)][is.nan(sbp), sbp := NA]
pr[, dbp := rowMeans(cbind(d2, d3), na.rm = TRUE)][is.nan(dbp), dbp := NA]
pr[, bio_date := make_date(ok(hv807y, 2019, 2021), ok(hv807m, 1, 12), ok(hv807d, 1, 31))]
pr[, int_date := make_date(hv007, hv006, hv016)]
pr[, date := fcoalesce(bio_date, int_date)]
pr[, htn_measured := as.integer(sbp >= 140 | dbp >= 90)]
pr[, htn_report   := as.integer(sbp >= 140 | dbp >= 90 | bp_med == 1)]
pr[, ate_cat := cut(shb53, c(-Inf, 0, 2, 4, 8, Inf), labels = c("<1h","1-2h","3-4h","5-8h",">8h"))]

d <- pr[hv105 >= 18 & preg == 0 & !is.na(sbp) & !is.na(dbp)]
d[, age_grp := cut(hv105, c(17, 29, 39, 49, 59, 69, Inf))]
d[, sub_bmi := (hv104 == 2 & hv105 <= 49) | (hv104 == 1 & hv105 <= 54)]
cat("Analytic adults (all dates):", nrow(d), "\n")

# Exclude 17 adults with measurement dates during the COVID suspension (likely date errors)
d <- d[!format(date, "%Y_%m") %in% c("2020_04", "2020_10")]
cat("Analytic adults (after excluding Apr & Oct 2020):", nrow(d), "\n")
saveRDS(d, file.path(OUT_DIR, "analytic_bp.rds"))
rm(pr); invisible(gc())

## =========================== 2. GPS ==========================================
ge <- st_read(GE_FILE, quiet = TRUE)
ge <- ge[!(ge$LATNUM == 0 & ge$LONGNUM == 0), c("DHSCLUST", "URBAN_RURA")]
cat("Clusters with GPS:", nrow(ge), "\n")

fill_coast <- function(r) cover(focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE),
                                focal(r, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE))

## ================= 3. CHECK DAILY CLIMATE FILES ==============================
WINDOW <- sort(sub("^t2m_(\\d{4}_\\d{2})\\.nc$", "\\1",
                   list.files(CLIM_DIR, pattern = "^t2m_\\d{4}_\\d{2}\\.nc$")))
WINDOW <- WINDOW[WINDOW >= "2019_05"]
files  <- file.path(CLIM_DIR, paste0("t2m_", WINDOW, ".nc"))
ref    <- rast(files[1])

qc <- rbindlist(lapply(files, function(f) {
  r   <- rast(f)
  dts <- as.Date(time(r))
  data.table(file = basename(f),
             same_grid = compareGeom(r, ref, crs = FALSE, stopOnError = FALSE),
             days = nlyr(r),
             days_expected = days_in_month(dts[1]),
             first = min(dts), last = max(dts),
             mean_C = round(global(mean(r), "mean", na.rm = TRUE)[[1]] - 273.15, 1))
}))
print(qc)
fwrite(qc, file.path(OUT_DIR, "climate_file_check.csv"))
if (any(!qc$same_grid) || any(qc$days != qc$days_expected))
  stop("Climate file problem: see table above")

## ============== 4. CLIMATOLOGICAL NORMALS (1991-2018) ========================
rmn <- rast(file.path(CLIM_DIR, "t2m_monthly_1991_2018.nc"))

# ---- diagnose ----
cat("Layers:", nlyr(rmn), "\n")
cat("Variables:", paste(varnames(rmn), collapse = ", "), "\n")
tt <- time(rmn)
cat("Time range:", format(min(tt, na.rm = TRUE)), "to", format(max(tt, na.rm = TRUE)),
    "| missing times:", sum(is.na(tt)), "\n")
print(head(names(rmn), 6))

# ---- keep one layer per month, 1991-2018, 2 m temperature only ----
if (length(varnames(rmn)) > 1) rmn <- rmn[[grep("t2m", names(rmn))]]
tt   <- time(rmn)
ym   <- format(tt, "%Y-%m")
yr   <- as.integer(format(tt, "%Y"))
keep <- which(!is.na(tt) & yr >= 1991 & yr <= 2018 & !duplicated(ym))
rmn  <- rmn[[keep]]
cat("Layers kept:", nlyr(rmn), "\n")
cat("Mean temperature of first layer (°C):",
    round(global(rmn[[1]], "mean", na.rm = TRUE)[[1]] - 273.15, 1), "\n")
stopifnot(nlyr(rmn) == 336)

chk <- vect(data.frame(city = c("Delhi", "Chennai"),
                       lon  = c(77.21, 80.27), lat = c(28.61, 13.08)),
            geom = c("lon", "lat"), crs = "EPSG:4326")
v   <- as.matrix(terra::extract(fill_coast(rmn), project(chk, crs(rmn)), ID = FALSE)) - 273.15
mon <- as.integer(format(time(rmn), "%m"))
print(data.frame(city = chk$city,                       # (repository fix: round numeric columns only)
                 round(data.frame(January = rowMeans(v[, mon == 1]),
                                  May     = rowMeans(v[, mon == 5]),
                                  July    = rowMeans(v[, mon == 7])), 1)))

if (FILL_NORMALS) rmn <- fill_coast(rmn)
pts_n <- project(vect(ge), crs(rmn))
mon   <- as.integer(format(time(rmn), "%m"))
vm    <- as.matrix(terra::extract(rmn, pts_n, ID = FALSE)) - 273.15
clim12 <- sapply(1:12, function(m) rowMeans(vm[, mon == m, drop = FALSE], na.rm = TRUE))

H   <- function(doy) cbind(1, sin(2*pi*doy/365.25), cos(2*pi*doy/365.25),
                           sin(4*pi*doy/365.25), cos(4*pi*doy/365.25))
mid <- c(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
B   <- solve(crossprod(H(mid)), t(H(mid)) %*% t(clim12))
norm_daily  <- H(1:366) %*% B
annual_mean <- data.table(DHSCLUST = ge$DHSCLUST, t_annual = rowMeans(clim12))
nd <- list(norm_daily = norm_daily, annual_mean = annual_mean, ids = ge$DHSCLUST)
saveRDS(nd, file.path(OUT_DIR, "normals.rds"))
cat("Clusters without normals:", sum(is.na(annual_mean$t_annual)), "\n")

## ============ 5. DAILY EXPOSURE AT CLUSTERS (full period) ====================
# Keep only cluster-days needed: measurement day, 1-2 days before, 7 days after
keys <- unique(d[, .(DHSCLUST = hv001, date)])
keys <- unique(rbindlist(lapply(c(-2, -1, 0, 7), function(k) keys[, .(DHSCLUST, date = date + k)])))

pts <- project(vect(ge), crs(ref))
daily <- rbindlist(lapply(files, function(f) {
  r   <- rast(f)
  dts <- as.Date(time(r))
  v   <- as.matrix(terra::extract(fill_coast(r), pts, ID = FALSE)) - 273.15
  dt  <- data.table(DHSCLUST = rep(ge$DHSCLUST, ncol(v)),
                    date = rep(dts, each = nrow(v)), t = as.vector(v))
  dt[keys, on = .(DHSCLUST, date), nomatch = NULL]
}))
raw <- terra::extract(ref[[1]], pts, ID = FALSE)[, 1]
daily[, filled := DHSCLUST %in% ge$DHSCLUST[is.na(raw)]]
daily[, anom := t - nd$norm_daily[cbind(yday(date), match(DHSCLUST, nd$ids))]]

for (k in 1:2) {
  tmp <- daily[, .(DHSCLUST, date = date + k, v = t)]; setnames(tmp, "v", paste0("t_lag", k))
  daily <- tmp[daily, on = .(DHSCLUST, date)]
}
tmp <- daily[, .(DHSCLUST, date = date - 7, t_lead7 = t)]
daily <- tmp[daily, on = .(DHSCLUST, date)]
saveRDS(daily, file.path(OUT_DIR, "daily_exposure.rds"))
cat("Clusters still missing temperature:", daily[is.na(t), uniqueN(DHSCLUST)], "\n")

## ======================= 6. ANALYSIS DATASET =================================
a <- merge(d, daily[, .(hv001 = DHSCLUST, date, t0 = t, anom0 = anom,
                        t_lag1, t_lag2, t_lead7, filled)], by = c("hv001", "date"))
a <- a[!is.na(t0)]
a <- merge(a, nd$annual_mean, by.x = "hv001", by.y = "DHSCLUST", all.x = TRUE)
a[, `:=`(ym = format(date, "%Y-%m"), dow = wday(date),
         t_mean02 = rowMeans(cbind(t0, t_lag1, t_lag2)),
         tbin = relevel(cut(t0, c(-Inf, 10, 15, 20, 25, 30, Inf),
                            labels = c("<10","10-15","15-20","20-25","25-30",">30")), ref = "20-25"),
         season = factor(fcase(month(date) %in% c(12, 1, 2), "winter",
                               month(date) %in% 3:5,         "pre-monsoon",
                               month(date) %in% 6:9,         "monsoon",
                               default = "post-monsoon"),
                         levels = c("winter", "pre-monsoon", "monsoon", "post-monsoon")))]
saveRDS(a, file.path(OUT_DIR, "analysis_dataset.rds"))
cat("Analysis N:", nrow(a), "| clusters:", uniqueN(a$hv001), "| districts:", uniqueN(a$shdist),
    "| states:", uniqueN(a$hv024), "\n")
print(a[, .N, by = season])

###############################################################################
# SECTIONS 7-10 (memory-safe): models, tables, figures
# Run after Session -> Restart R. Uses output_full/analysis_dataset.rds
###############################################################################
library(data.table); library(fixest); library(ggplot2)

OUT_DIR <- "output_full"
state_lab <- c("1"="Jammu & Kashmir","2"="Himachal Pradesh","3"="Punjab","4"="Chandigarh",
               "5"="Uttarakhand","6"="Haryana","7"="Delhi","8"="Rajasthan","9"="Uttar Pradesh","10"="Bihar",
               "11"="Sikkim","12"="Arunachal Pradesh","13"="Nagaland","14"="Manipur","15"="Mizoram",
               "16"="Tripura","17"="Meghalaya","18"="Assam","19"="West Bengal","20"="Jharkhand","21"="Odisha",
               "22"="Chhattisgarh","23"="Madhya Pradesh","24"="Gujarat","25"="Dadra & Nagar Haveli and Daman & Diu",
               "27"="Maharashtra","28"="Andhra Pradesh","29"="Karnataka","30"="Goa","31"="Lakshadweep",
               "32"="Kerala","33"="Tamil Nadu","34"="Puducherry","35"="Andaman & Nicobar Islands",
               "36"="Telangana","37"="Ladakh")

## ---------------- load analysis data and keep only needed columns -----------
a <- readRDS(file.path(OUT_DIR, "analysis_dataset.rds"))
keep_cols <- c("hv001","hv024","shdist","shphase","date","ym","dow","season","wt",
               "hv104","age_grp","hv270","hv025","shb14a","shb14b","shb14c","sh25","sh26",
               "sbp","dbp","htn_measured","htn_report","glu","glu_med","ate_cat","bp_med",
               "height","bmi","sub_bmi","hv040","filled","fan","cooler",
               "t0","anom0","tbin","t_lead7","t_mean02","t_annual")
a <- a[, ..keep_cols]
invisible(gc())
cat("Analysis N:", nrow(a), "| clusters:", uniqueN(a$hv001), "| districts:", uniqueN(a$shdist),
    "| states:", uniqueN(a$hv024), "\n")
print(a[, .N, by = season])

## ============================ 7. MODELS ======================================
ctrl     <- paste("+ age_grp*hv104 + i(hv270) + hv025 + i(shb14a) + i(shb14b) + i(shb14c)",
                  "+ i(sh25) + i(sh26)")
rhs      <- paste(ctrl, "| shdist + ym + dow")
rhs_dxym <- paste(ctrl, "| shdist^ym + dow")

fit <- function(lhs, x, sub = NULL, r = rhs, cl = ~hv001) {
  f <- as.formula(paste(lhs, "~", x, r))
  mod <- if (is.null(sub)) {
    feols(f, a, weights = ~wt, cluster = cl, lean = TRUE, mem.clean = TRUE)
  } else {
    feols(f, a, weights = ~wt, cluster = cl, subset = sub, lean = TRUE, mem.clean = TRUE)
  }
  invisible(gc())
  mod
}

specs <- list(
  sbp      = list(lhs = "sbp",          x = "t0"),
  dbp      = list(lhs = "dbp",          x = "t0"),
  htn      = list(lhs = "htn_measured", x = "t0"),
  glu      = list(lhs = "log(glu)",     x = "t0 + i(ate_cat)", sub = ~ glu_med == 0),
  anom     = list(lhs = "sbp",          x = "anom0"),
  bins     = list(lhs = "sbp",          x = "i(tbin, ref = '20-25')"),
  lead     = list(lhs = "sbp",          x = "t_lead7"),
  joint    = list(lhs = "sbp",          x = "t0 + t_lead7"),
  height   = list(lhs = "height",       x = "t0", sub = ~ sub_bmi == TRUE),
  nomeds   = list(lhs = "sbp",          x = "t0", sub = ~ bp_med == 0),
  noalt    = list(lhs = "sbp",          x = "t0", sub = ~ hv040 <= 2000 & filled == FALSE),
  mean02   = list(lhs = "sbp",          x = "t_mean02"),
  bmi      = list(lhs = "sbp",          x = "t0 + bmi", sub = ~ sub_bmi == TRUE),
  sex      = list(lhs = "sbp",          x = "i(hv104, t0)"),
  age      = list(lhs = "sbp",          x = "i(age_grp, t0)"),
  cooling  = list(lhs = "sbp",          x = "t0 * fan + t0 * cooler"),
  season   = list(lhs = "sbp",          x = "i(season, t0)"),
  dxym     = list(lhs = "sbp",          x = "t0", r = rhs_dxym),
  dxym_dbp = list(lhs = "dbp",          x = "t0", r = rhs_dxym),
  sbp_dcl  = list(lhs = "sbp",          x = "t0", cl = ~shdist)
)

# Fit one model at a time; save after each so a crash never loses finished models
mfile <- file.path(OUT_DIR, "models.rds")
m <- if (file.exists(mfile)) readRDS(mfile) else list()
for (nm in setdiff(names(specs), names(m))) {
  m[[nm]] <- do.call(fit, specs[[nm]])
  saveRDS(m, mfile)
  cat(nm, "done | N =", nobs(m[[nm]]), "\n")
}

etable(m$sbp, m$dbp, m$htn, m$glu, m$anom, m$dxym, keep = c("t0", "anom0"))

## ============================ 8. TABLES ======================================
w   <- function(x) round(100 * weighted.mean(x, a$wt, na.rm = TRUE), 1)
wm  <- function(x) round(weighted.mean(x, a$wt, na.rm = TRUE), 1)
wsd <- function(x) round(sqrt(Hmisc::wtd.var(x, a$wt, na.rm = TRUE)), 1)

# Table 1: sample characteristics
tab1 <- rbind(
  data.table(item = "N; clusters; districts; states",
             value = paste(nrow(a), uniqueN(a$hv001), uniqueN(a$shdist), uniqueN(a$hv024), sep = "; ")),
  a[, .(item = paste("Age", age_grp), value = round(100 * sum(wt) / a[, sum(wt)], 1)), by = age_grp
  ][order(age_grp)][, .(item, value = as.character(value))],
  a[, .(item = paste("Season:", season), value = round(100 * sum(wt) / a[, sum(wt)], 1)), by = season
  ][order(season)][, .(item, value = as.character(value))],
  data.table(item = c("Women %", "Urban %", "Tobacco %", "Alcohol %", "SBP mean (SD)", "DBP mean (SD)",
                      "Measured HTN %", "HTN incl. meds %", "Glucose mean mg/dL",
                      "Temp mean (IQR; P5-P95)"),
             value = c(w(a$hv104 == 2), w(a$hv025 == 1), w(a$sh25 == 1), w(a$sh26 == 1),
                       paste0(wm(a$sbp), " (", wsd(a$sbp), ")"), paste0(wm(a$dbp), " (", wsd(a$dbp), ")"),
                       w(a$htn_measured), w(a$htn_report), wm(a$glu),
                       sprintf("%.1f (%.1f-%.1f; %.1f-%.1f)", mean(a$t0),
                               quantile(a$t0, .25), quantile(a$t0, .75),
                               quantile(a$t0, .05), quantile(a$t0, .95)))))
tab1 <- rbind(tab1, a[, .(item = paste("Wealth Q", hv270),
                          value = as.character(round(100 * sum(wt) / a[, sum(wt)], 1))), by = hv270
][order(hv270)][, .(item, value)])
fwrite(tab1, file.path(OUT_DIR, "table1.csv")); print(tab1)

# Table 2: estimates per 1 deg C with 95% CI
ci_row <- function(label, mod, term, scale = 1) {
  ct  <- coeftable(mod)
  idx <- if (term %in% rownames(ct)) term else grep(term, rownames(ct), fixed = TRUE)[1]
  b <- ct[idx, 1]; se <- ct[idx, 2]
  data.table(model = label, n = nobs(mod), est = round(scale * b, 3),
             lo = round(scale * (b - 1.96 * se), 3), hi = round(scale * (b + 1.96 * se), 3))
}
age_terms    <- grep("t0", rownames(coeftable(m$age)), value = TRUE)
season_terms <- grep("season::", rownames(coeftable(m$season)), value = TRUE)

tab2 <- rbindlist(list(
  ci_row("SBP, mmHg", m$sbp, "t0"),
  ci_row("DBP, mmHg", m$dbp, "t0"),
  ci_row("Measured hypertension, pp", m$htn, "t0", 100),
  ci_row("Random glucose, %", m$glu, "t0", 100),
  ci_row("SBP: anomaly", m$anom, "anom0"),
  ci_row("SBP: mean days 0-2", m$mean02, "t_mean02"),
  ci_row("SBP: excl. BP meds", m$nomeds, "t0"),
  ci_row("SBP: excl. >2000 m & coast-filled", m$noalt, "t0"),
  ci_row("SBP: + BMI", m$bmi, "t0"),
  ci_row("SBP: SE clustered by district", m$sbp_dcl, "t0"),
  ci_row("SBP: district x month FE", m$dxym, "t0"),
  ci_row("DBP: district x month FE", m$dxym_dbp, "t0"),
  ci_row("Placebo alone: temp 7 d later", m$lead, "t_lead7"),
  ci_row("Joint: same-day temp", m$joint, "t0"),
  ci_row("Joint: temp 7 d later", m$joint, "t_lead7"),
  ci_row("Negative control: height, cm", m$height, "t0"),
  ci_row("SBP: men", m$sex, "hv104::1"),
  ci_row("SBP: women", m$sex, "hv104::2"),
  rbindlist(lapply(age_terms,    function(tt) ci_row(paste("SBP:", tt), m$age, tt))),
  rbindlist(lapply(season_terms, function(tt) ci_row(paste("SBP:", tt), m$season, tt))),
  ci_row("SBP: t0 x fan (interaction)", m$cooling, "t0:fan")
))
fwrite(tab2, file.path(OUT_DIR, "table2.csv")); print(tab2, nrows = 60)

# Table 3: temperature-standardised measured hypertension by state
b_s <- coef(m$sbp)["t0"]; b_d <- coef(m$dbp)["t0"]
s <- a[!is.na(t_annual)]
s[, `:=`(sbp_adj = sbp - b_s * (t0 - t_annual), dbp_adj = dbp - b_d * (t0 - t_annual))]
s[, htn_adj := as.integer(sbp_adj >= 140 | dbp_adj >= 90)]

tab3 <- s[, .(n = .N, mean_t0 = round(mean(t0), 1), mean_annual = round(mean(t_annual), 1),
              obs = 100 * weighted.mean(htn_measured, wt),
              adj = 100 * weighted.mean(htn_adj, wt)), by = hv024][n >= 500]
tab3[, `:=`(state = state_lab[as.character(hv024)], change_pp = adj - obs,
            rank_obs = frank(-obs), rank_adj = frank(-adj))]
tab3 <- tab3[order(rank_obs), .(state, n, mean_t0, mean_annual, obs = round(obs, 1), adj = round(adj, 1),
                                change_pp = round(change_pp, 1), rank_obs, rank_adj)]
pooled <- s[, .(obs = 100 * weighted.mean(htn_measured, wt), adj = 100 * weighted.mean(htn_adj, wt))]
fwrite(tab3, file.path(OUT_DIR, "table3_state_standardised.csv"))
print(tab3, nrows = 50); print(round(pooled, 1))
cat("States changing rank:", tab3[rank_obs != rank_adj, .N], "of", nrow(tab3), "\n")
rm(s); invisible(gc())

## ============================ 9. FIGURES =====================================
# Figure 1: fieldwork calendar
cal <- a[, .N, by = .(month = as.Date(format(date, "%Y-%m-01")), phase = factor(shphase))]
f1 <- ggplot(cal, aes(month, N / 1000, fill = phase)) +
  geom_col() + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Month of biomarker measurement", y = "Adults measured (thousands)", fill = "Phase") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "fig1_calendar.png"), f1, width = 8, height = 4, dpi = 300)

# Figure 2: SBP by temperature category (ref 20-25 C)
ct <- as.data.table(coeftable(m$bins), keep.rownames = "term")[grepl("^tbin::", term)]
ct[, bin := sub("tbin::", "", term)]
f2d <- rbind(ct[, .(bin, est = Estimate, lo = Estimate - 1.96 * `Std. Error`,
                    hi = Estimate + 1.96 * `Std. Error`)],
             data.table(bin = "20-25", est = 0, lo = 0, hi = 0))
f2d[, bin := factor(bin, levels = c("<10","10-15","15-20","20-25","25-30",">30"))]
f2 <- ggplot(f2d, aes(bin, est)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  labs(x = "Measurement-day temperature (°C)",
       y = "Difference in systolic BP (mmHg)\nvs 20–25 °C (95% CI)") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "fig2_temperature_bins.png"), f2, width = 7, height = 4.5, dpi = 300)

# Figure 3: SBP-temperature slope by season
sd3 <- as.data.table(coeftable(m$season), keep.rownames = "term")[grepl("season::", term)]
sd3[, season := factor(sub(".*season::([a-z-]+).*", "\\1", term),
                       levels = c("winter", "pre-monsoon", "monsoon", "post-monsoon"))]
f3 <- ggplot(sd3, aes(season, Estimate)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = Estimate - 1.96 * `Std. Error`, ymax = Estimate + 1.96 * `Std. Error`)) +
  labs(x = "Season of measurement", y = "Change in systolic BP (mmHg)\nper 1 °C (95% CI)") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "fig3_season_slopes.png"), f3, width = 7, height = 4.5, dpi = 300)

## ========================= 10. RECORD SESSION ================================
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("Done. Outputs in", OUT_DIR, "\n")

## ---- A. Do the season slopes really differ? ----
b  <- coef(m$season); V <- vcov(m$season)
nm <- grep("season::", names(b), value = TRUE)
nm <- nm[order(match(sub(".*season::([a-z-]+).*", "\\1", nm),
                     c("winter", "pre-monsoon", "monsoon", "post-monsoon")))]
print(nm)

pair_test <- function(t1, t2) {
  d  <- b[t1] - b[t2]
  se <- sqrt(V[t1, t1] + V[t2, t2] - 2 * V[t1, t2])
  data.table(comparison = paste(sub(".*season::([a-z-]+).*", "\\1", t1), "vs",
                                sub(".*season::([a-z-]+).*", "\\1", t2)),
             diff = round(d, 3), lo = round(d - 1.96 * se, 3), hi = round(d + 1.96 * se, 3),
             p = signif(2 * pnorm(-abs(d / se)), 3))
}
season_tests <- rbind(pair_test(nm[3], nm[1]), pair_test(nm[2], nm[1]), pair_test(nm[4], nm[1]))
print(season_tests)

# Joint test: are all four season slopes equal?
L  <- rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1))
Lb <- L %*% b[nm]
W  <- as.numeric(t(Lb) %*% solve(L %*% V[nm, nm] %*% t(L)) %*% Lb)
cat("Joint Wald test of equal season slopes: chi2 =", round(W, 1), "on 3 df, p =",
    signif(pchisq(W, 3, lower.tail = FALSE), 3), "\n")
fwrite(season_tests, file.path(OUT_DIR, "season_slope_tests.csv"))

## ---- B. Table 3 sensitivity: clusters at or below 2000 m only ----
b_s <- coef(m$sbp)["t0"]; b_d <- coef(m$dbp)["t0"]
s <- a[!is.na(t_annual) & hv040 <= 2000]
s[, `:=`(sbp_adj = sbp - b_s * (t0 - t_annual), dbp_adj = dbp - b_d * (t0 - t_annual))]
s[, htn_adj := as.integer(sbp_adj >= 140 | dbp_adj >= 90)]

tab3b <- s[, .(n = .N, mean_t0 = round(mean(t0), 1), mean_annual = round(mean(t_annual), 1),
               obs = 100 * weighted.mean(htn_measured, wt),
               adj = 100 * weighted.mean(htn_adj, wt)), by = hv024][n >= 500]
tab3b[, `:=`(state = state_lab[as.character(hv024)], change_pp = adj - obs,
             rank_obs = frank(-obs), rank_adj = frank(-adj))]
tab3b <- tab3b[order(rank_obs), .(state, n, mean_t0, mean_annual, obs = round(obs, 1),
                                  adj = round(adj, 1), change_pp = round(change_pp, 1),
                                  rank_obs, rank_adj)]
pooled_b <- s[, .(obs = 100 * weighted.mean(htn_measured, wt), adj = 100 * weighted.mean(htn_adj, wt))]
fwrite(tab3b, file.path(OUT_DIR, "table3_state_standardised_le2000m.csv"))
print(tab3b, nrows = 50); print(round(pooled_b, 1))
cat("States moving 3+ ranks:", tab3b[abs(rank_obs - rank_adj) >= 3, .N], "of", nrow(tab3b), "\n")

## ---- C. Pooled line for the main Table 3 (all clusters) ----
s_all <- a[!is.na(t_annual)]
s_all[, htn_adj := as.integer((sbp - b_s * (t0 - t_annual)) >= 140 | (dbp - b_d * (t0 - t_annual)) >= 90)]
print(round(s_all[, .(obs = 100 * weighted.mean(htn_measured, wt),
                      adj = 100 * weighted.mean(htn_adj, wt))], 1))
rm(s, s_all); invisible(gc())


## 1. Exact values for the temperature-bin results (replace the [verify] numbers)
round(coeftable(m$bins)[, 1:2], 3)

## 2. Confirm Lakshadweep (code 31) is absent, and get the temperature range
cat("Lakshadweep in sample:", 31 %in% a$hv024, "\n")
print(round(range(a$t0), 1))

## 3. Table 3 direction check: weighted mean shift and who crosses the threshold
s <- a[!is.na(t_annual)]
s[, shift := t0 - t_annual]
s[, htn_adj := as.integer((sbp - b_s * shift) >= 140 | (dbp - b_d * shift) >= 90)]
chk <- s[, .(n = .N,
             w_shift  = round(weighted.mean(shift, wt), 2),
             obs      = round(100 * weighted.mean(htn_measured, wt), 2),
             adj      = round(100 * weighted.mean(htn_adj, wt), 2),
             moved_up = round(100 * weighted.mean(htn_adj == 1 & htn_measured == 0, wt), 2),
             moved_dn = round(100 * weighted.mean(htn_adj == 0 & htn_measured == 1, wt), 2)),
         by = hv024]
chk[, state := state_lab[as.character(hv024)]]
print(chk[order(w_shift), .(state, n, w_shift, obs, adj, moved_up, moved_dn)], nrows = 50)
rm(s); invisible(gc())


## ---- Diagnostic: values are on a 0.5 mm Hg grid and heaped at the thresholds ----
cat("SBP on 0.5 grid:", round(100 * a[, mean(sbp %% 0.5 == 0)], 2), "% | DBP on 0.5 grid:",
    round(100 * a[, mean(dbp %% 0.5 == 0)], 2), "%\n")
cat("Weighted % with SBP exactly 140.0:", round(100 * a[, weighted.mean(sbp == 140, wt)], 2),
    "| exactly 139.5:", round(100 * a[, weighted.mean(sbp == 139.5, wt)], 2), "\n")
cat("Weighted % with DBP exactly 90.0:", round(100 * a[, weighted.mean(dbp == 90, wt)], 2),
    "| exactly 89.5:", round(100 * a[, weighted.mean(dbp == 89.5, wt)], 2), "\n")

## ---- Corrected standardisation ----
b_s  <- coef(m$sbp)["t0"]; b_d <- coef(m$dbp)["t0"]
snap <- function(x) round(2 * x) / 2          # back onto the 0.5 mm Hg measurement grid

std_table <- function(dat) {
  s <- dat[!is.na(t_annual)]
  s[, shift := t0 - t_annual]
  s[, `:=`(sbp_adj = snap(sbp - b_s * shift), dbp_adj = snap(dbp - b_d * shift))]
  s[, htn_adj := as.integer(sbp_adj >= 140 | dbp_adj >= 90)]
  tab <- s[, .(n = .N, mean_t0 = round(mean(t0), 1), mean_annual = round(mean(t_annual), 1),
               w_shift  = round(weighted.mean(shift, wt), 2),
               obs      = 100 * weighted.mean(htn_measured, wt),
               adj      = 100 * weighted.mean(htn_adj, wt),
               moved_up = round(100 * weighted.mean(htn_adj == 1 & htn_measured == 0, wt), 2),
               moved_dn = round(100 * weighted.mean(htn_adj == 0 & htn_measured == 1, wt), 2)),
           by = hv024][n >= 500]
  tab[, `:=`(state = state_lab[as.character(hv024)], change_pp = adj - obs,
             rank_obs = frank(-obs), rank_adj = frank(-adj))]
  tab <- tab[order(rank_obs), .(state, n, mean_t0, mean_annual, w_shift,
                                obs = round(obs, 1), adj = round(adj, 1),
                                change_pp = round(change_pp, 1), rank_obs, rank_adj,
                                moved_up, moved_dn)]
  pooled <- s[, .(obs = round(100 * weighted.mean(htn_measured, wt), 1),
                  adj = round(100 * weighted.mean(htn_adj, wt), 1))]
  list(tab = tab, pooled = pooled)
}

main <- std_table(a)
low  <- std_table(a[hv040 <= 2000])

print(main$tab, nrows = 50); print(main$pooled)
cat("States changing rank:", main$tab[rank_obs != rank_adj, .N],
    "| moving 3+ ranks:", main$tab[abs(rank_obs - rank_adj) >= 3, .N],
    "of", nrow(main$tab), "\n")
cat("States where prevalence fell:", main$tab[change_pp < 0, .N],
    "| rose:", main$tab[change_pp > 0, .N], "\n")
print(low$pooled)

fwrite(main$tab, file.path(OUT_DIR, "table3_state_standardised_v2.csv"))
fwrite(low$tab,  file.path(OUT_DIR, "table3_state_standardised_le2000m_v2.csv"))

print(main$pooled)

library(data.table); library(fixest)
OUT_DIR <- "output_full"
a <- readRDS(file.path(OUT_DIR, "analysis_dataset.rds"))
a <- a[, .(hv001, hv024, shdist, ym, dow, wt, hv104, age_grp, hv270, hv025,
           shb14a, shb14b, shb14c, sh25, sh26, sbp, dbp, t0, bp_hour, hv018)]
invisible(gc())

# Hour of first BP reading, grouped; missing time kept as its own category
a[, hour_cat := as.character(cut(bp_hour, c(6, 9, 11, 13, 15, 17, 21), include.lowest = TRUE))]
a[is.na(hour_cat), hour_cat := "missing"]
print(a[, .N, by = hour_cat][order(hour_cat)])

ctrl <- "+ age_grp*hv104 + i(hv270) + hv025 + i(shb14a) + i(shb14b) + i(shb14c) + i(sh25) + i(sh26)"

m_hour <- feols(as.formula(paste("sbp ~ t0", ctrl, "+ i(hour_cat) | shdist + ym + dow")),
                a, weights = ~wt, cluster = ~hv001, lean = TRUE, mem.clean = TRUE)
invisible(gc())
m_int  <- feols(as.formula(paste("sbp ~ t0", ctrl, "| shdist + ym + dow + hv024^hv018")),
                a, weights = ~wt, cluster = ~hv001, lean = TRUE, mem.clean = TRUE)

etable(m_hour, m_int, keep = "t0")
for (mm in list(m_hour, m_int)) {
  ct <- coeftable(mm)["t0", ]
  cat(sprintf("n = %d | t0 = %.3f (%.3f to %.3f)\n", nobs(mm),
              ct[1], ct[1] - 1.96 * ct[2], ct[1] + 1.96 * ct[2]))
}
