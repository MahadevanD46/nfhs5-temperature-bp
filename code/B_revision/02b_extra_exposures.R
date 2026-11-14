# =============================================================================
# 02b_extra_exposures.R — new exposures and variables for the revision
#   * interview-date temperature                               (comment 10)
#   * coastal infill method per cluster: own cell / 3x3 / 5x5   (comment 13)
#   * displacement-weighted temperature + cell-change probability (comment 12)
#   * dewpoint, relative humidity, precipitation                (comments 15, 19)
#   * ERA5-Land cell id, coordinates for spatial SEs            (comment 16)
#   * gap between interview and biomarker dates, date source    (comment 10)
# Output: output_revision/analysis_dataset_v2.rds  (RESTRICTED: never share)
#         output_revision/coastal_infill_by_state.csv          -> Suppl. Table S7
#         output_revision/displacement_summary.csv             -> Suppl. Table S10 note
# =============================================================================
source("code/B_revision/00_config.R")
suppressPackageStartupMessages({ library(terra); library(sf); library(Matrix); library(lubridate) })

a <- load_a(extra = FALSE)
keep_cols <- c("hv001","hv002","hvidx","hv021","hv022","hv024","shdist","shphase","date",
  "bio_date","int_date","ym","dow","season","wt","hv104","hv105","age_grp","hv270","hv025",
  "shb14a","shb14b","shb14c","sh25","sh26","sbp","dbp","htn_measured","htn_report","bp_med",
  "glu","glu_med","ate_cat","height","bmi","sub_bmi","hv040","filled","fan","cooler",
  "bp_hour","hv018","t0","anom0","tbin","t_lead7","t_mean02","t_lag1","t_lag2","t_annual")
a <- a[, intersect(keep_cols, names(a)), with = FALSE]; invisible(gc())

## ---- clusters and grid ------------------------------------------------------
ge <- st_read(GE_FILE, quiet = TRUE)
ge <- ge[!(ge$LATNUM == 0 & ge$LONGNUM == 0), c("DHSCLUST", "URBAN_RURA", "LATNUM", "LONGNUM")]
gdt <- as.data.table(st_drop_geometry(ge))[DHSCLUST %in% unique(a$hv001)]
files <- sort(list.files(CLIM_DIR, pattern = "^t2m_\\d{4}_\\d{2}\\.nc$", full.names = TRUE))
files <- files[sub("^t2m_(\\d{4}_\\d{2})\\.nc$", "\\1", basename(files)) >= "2019_05"]
ref <- rast(files[1])
fill_coast <- function(r) cover(focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE),
                                focal(r, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE))
xy <- cbind(gdt$LONGNUM, gdt$LATNUM)
gdt[, cell := cellFromXY(ref, xy)]

## ---- (13) which infill step was used ---------------------------------------
r1 <- ref[[1]]
v_raw <- terra::extract(r1, xy)[, 1]
v_3   <- terra::extract(focal(r1, 3, "mean", na.policy = "only", na.rm = TRUE), xy)[, 1]
v_5   <- terra::extract(focal(r1, 5, "mean", na.policy = "only", na.rm = TRUE), xy)[, 1]
gdt[, coast_method := fcase(!is.na(v_raw), "own cell", !is.na(v_3), "3x3",
                            !is.na(v_5), "5x5", default = "unfilled")]

## ---- (12) displacement simulation ------------------------------------------
# DHS rule: random direction, distance uniform on [0, R]; R = 2 km urban, 5 km rural
# (1% of rural clusters up to 10 km; ignored, so estimates are slightly optimistic).
K <- SET$n_disp_draws
sim <- gdt[, .(k = seq_len(K)), by = .(DHSCLUST, LATNUM, LONGNUM, URBAN_RURA, cell)]
sim[, R := fifelse(URBAN_RURA == "U", 2, 5)]
sim[, `:=`(ang = runif(.N, 0, 2 * pi), dist = runif(.N, 0, R))]
sim[, cell2 := cellFromXY(ref, cbind(LONGNUM + dist * cos(ang) / (111.32 * cos(LATNUM * pi / 180)),
                                    LATNUM + dist * sin(ang) / 111.32))]
disp <- sim[, .(p_same_cell = mean(cell2 == cell, na.rm = TRUE),
                n_cells = uniqueN(cell2)), by = DHSCLUST]
gdt <- merge(gdt, disp, by = "DHSCLUST")
trip <- sim[!is.na(cell2), .(w = .N), by = .(DHSCLUST, cell = cell2)]
cells_needed <- sort(unique(c(trip$cell, gdt$cell)))
gdt[, row := .I]
W <- sparseMatrix(i = gdt$row[match(trip$DHSCLUST, gdt$DHSCLUST)],
                  j = match(trip$cell, cells_needed), x = trip$w,
                  dims = c(nrow(gdt), length(cells_needed)))
rm(sim); invisible(gc())

## ---- loop over months: interview-date and displacement-weighted temperature ---
keys_m <- unique(a[, .(hv001, date)])
keys_i <- unique(a[!is.na(int_date), .(hv001, date = int_date)])
own_j  <- match(gdt$cell, cells_needed)
res_m <- list(); res_i <- list()
for (f in files) {
  r <- fill_coast(rast(f)); dts <- as.Date(time(r))
  V <- as.matrix(terra::extract(r, cells_needed)) - 273.15           # cells x days
  ok <- !is.na(V)
  V0 <- V; V0[!ok] <- 0
  Td <- as.matrix(W %*% V0) / as.matrix(W %*% (ok * 1))              # renormalise over land
  To <- V[own_j, , drop = FALSE]                                       # = original t0 rule
  long <- function(M) data.table(hv001 = rep(gdt$DHSCLUST, ncol(M)),
                                 date = rep(dts, each = nrow(M)), v = as.vector(M))
  res_m[[f]] <- long(Td)[keys_m, on = .(hv001, date), nomatch = NULL]
  res_i[[f]] <- long(To)[keys_i, on = .(hv001, date), nomatch = NULL]
  log_msg(basename(f)); rm(V, V0, Td, To, r); invisible(gc())
}
t_disp <- rbindlist(res_m); setnames(t_disp, "v", "t_disp")
t_int  <- rbindlist(res_i); setnames(t_int, c("date", "v"), c("int_date", "t_int"))

## ---- dewpoint and precipitation at the measurement date ----------------------
extract_var <- function(pattern, varname, transform) {
  ff <- sort(list.files(CLIM_DIR, pattern = pattern, full.names = TRUE))
  if (!length(ff)) { log_msg("no files for", varname, "- skipped"); return(NULL) }
  pts <- cbind(gdt$LONGNUM, gdt$LATNUM)
  x <- rbindlist(lapply(ff, function(f) {
    r <- fill_coast(rast(f)); dts <- as.Date(time(r))
    v <- transform(as.matrix(terra::extract(r, pts)))
    data.table(hv001 = rep(gdt$DHSCLUST, ncol(v)), date = rep(dts, each = nrow(v)),
               v = as.vector(v))[keys_m, on = .(hv001, date), nomatch = NULL]
  }))
  setnames(x, "v", varname); x
}
d2m <- extract_var("^d2m_\\d{4}_\\d{2}\\.nc$", "d2m", function(x) x - 273.15)
tp  <- extract_var("^tp_\\d{4}_\\d{2}\\.nc$", "tp_mm", function(x) pmax(x * 1000, 0))

## ---- merge and derive -------------------------------------------------------
a <- merge(a, t_disp, by = c("hv001", "date"), all.x = TRUE)
a <- merge(a, t_int, by = c("hv001", "int_date"), all.x = TRUE)
if (!is.null(d2m)) a <- merge(a, d2m, by = c("hv001", "date"), all.x = TRUE)
if (!is.null(tp))  a <- merge(a, tp,  by = c("hv001", "date"), all.x = TRUE)
a <- merge(a, gdt[, .(hv001 = DHSCLUST, cell, coast_method, p_same_cell, n_cells,
                      lat = LATNUM, lon = LONGNUM)], by = "hv001", all.x = TRUE)
a[is.na(t_disp), t_disp := t0]       # all simulated locations at sea -> own (filled) cell
es <- function(T) 6.112 * exp(17.62 * T / (243.12 + T))              # Magnus (hPa)
if ("d2m" %in% names(a)) a[, rh := pmin(100, 100 * es(d2m) / es(t0))]
if ("tp_mm" %in% names(a)) a[, tp_cat := cut(tp_mm, c(-Inf, 0.1, 2.5, 7.5, 35.5, Inf),
                          labels = c("none", "very light", "light", "moderate", "heavy+"))]
a[, `:=`(gap = as.integer(date - int_date),
         same_day = as.integer(!is.na(bio_date) & bio_date == int_date),
         date_fallback = as.integer(is.na(bio_date)),
         t_doy = t0 - anom0,                                           # day-of-year normal
         shift = t0 - t_annual,
         region = region_of[as.character(hv024)],
         cell_date = paste(cell, date),
         hour_cat = as.character(cut(bp_hour, c(6, 9, 11, 13, 15, 17, 21), include.lowest = TRUE)))]
a[is.na(hour_cat), hour_cat := "missing"]
a[, gap_cat := cut(abs(gap), c(-1, 0, 1, 3, 7, Inf), labels = c("0", "1", "2-3", "4-7", ">7"))]

# season-mean normal for each cluster (alternative reference temperature, comment 22)
nd <- readRDS(file.path(IN_DIR, "normals.rds"))
doy_season <- function(doy) { m <- month(as.Date(doy - 1, origin = "2021-01-01"))
  fcase(m %in% c(12, 1, 2), "winter", m %in% 3:5, "pre-monsoon", m %in% 6:9, "monsoon",
        default = "post-monsoon") }
sea_of_doy <- doy_season(1:365)
snorm <- rbindlist(lapply(unique(sea_of_doy), function(s)
  data.table(hv001 = nd$ids, season = s,
             t_season_norm = colMeans(nd$norm_daily[which(sea_of_doy == s), , drop = FALSE]))))
a[, season_chr := as.character(season)]
a <- merge(a, snorm, by.x = c("hv001", "season_chr"), by.y = c("hv001", "season"), all.x = TRUE)
a[, season_chr := NULL]

saveRDS(a, file.path(OUT_DIR, "analysis_dataset_v2.rds"))
log_msg("analysis_dataset_v2.rds:", nrow(a), "rows")

## ---- summaries ---------------------------------------------------------------
gdt[, state := a$hv024[match(DHSCLUST, a$hv001)]]
pp <- a[, .N, by = .(hv001)]
gdt <- merge(gdt, pp, by.x = "DHSCLUST", by.y = "hv001", all.x = TRUE)
cs <- gdt[, .(clusters = .N, participants = sum(N, na.rm = TRUE)), by = .(coast_method, state)]
cs[, state := state_lab[as.character(state)]]
fwrite(cs[order(coast_method, -participants)], file.path(OUT_DIR, "coastal_infill_by_state.csv"))
tot <- gdt[, .(clusters = .N, participants = sum(N, na.rm = TRUE)), by = coast_method]
print(tot)
for (m in c("3x3", "5x5")) {
  k <- gsub("x", "X", m)
  reg(paste0("COAST_CL_", k), f0(tot[coast_method == m, sum(clusters)]))
  reg(paste0("COAST_PP_", k), f0(tot[coast_method == m, sum(participants)]))
}
reg("COAST_STATES", paste(cs[coast_method != "own cell", .(p = sum(participants)), by = state
                             ][order(-p)][1:min(.N, 6), state], collapse = ", "))
reg("NORMALS_MISSING", f0(a[is.na(t_annual), uniqueN(hv001)]))

ds <- gdt[, .(clusters = .N, pct_cross = 100 * mean(p_same_cell < 1),
              mean_p_other = 100 * mean(1 - p_same_cell)), by = .(urban = URBAN_RURA)]
ds <- rbind(ds, gdt[, .(urban = "All", clusters = .N, pct_cross = 100 * mean(p_same_cell < 1),
                        mean_p_other = 100 * mean(1 - p_same_cell))])
fwrite(ds, file.path(OUT_DIR, "displacement_summary.csv")); print(ds)
reg("DISP_PCT_CROSS", f1(ds[urban == "All", pct_cross]))
reg("DISP_PCT_CROSS_U", f1(ds[urban == "U", pct_cross]))
reg("DISP_PCT_CROSS_R", f1(ds[urban == "R", pct_cross]))
reg("DISP_P_OTHER", f1(ds[urban == "All", mean_p_other]))
reg("DISP_TDIFF_SD", f2(a[, sd(t_disp - t0, na.rm = TRUE)]))
reg("DISP_TCOR", formatC(a[, cor(t_disp, t0, use = "complete.obs")], format = "f", digits = 3))
log_msg("02b done")
