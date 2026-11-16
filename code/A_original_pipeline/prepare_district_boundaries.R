###############################################################################
# prepare_district_boundaries.R
# NFHS-5 district boundaries (DHS Spatial Data Repository, via the rdhs package)
# with one polygon per NFHS-5 district code (shdist).
#   Inputs : IAPR7EFL-derived analysis dataset (output_full/analysis_dataset.rds),
#            IAGE7AFL.shp (cluster GPS)
#   Output : sdr_boundaries/nfhs5_districts_fixed.gpkg  (field "shdist")
# Run once, after nfhs5_temperature_full_period.R and before the B_revision scripts.
# rdhs needs a DHS Program account (see ?rdhs::set_rdhs_config).
###############################################################################
library(rdhs); library(data.table); library(sf)
sf_use_s2(FALSE)
OUT_DIR <- "output_full"
dir.create("sdr_boundaries", showWarnings = FALSE)
raw_file <- "sdr_boundaries/nfhs5_districts.gpkg"

## 1. Download the district layer (first run only)
if (!file.exists(raw_file)) {
  sv  <- dhs_surveys(countryIds = "IA")
  sid <- sv$SurveyId[sv$SurveyYearLabel == "2019-21" & sv$SurveyType == "DHS"][1]
  bnd <- download_boundaries(surveyId = sid, method = "sf")
  dist <- bnd[[which.min(abs(sapply(bnd, nrow) - 707))]]          # the district layer (about 707 polygons)
  st_write(dist, raw_file, quiet = TRUE)
}
dist <- st_make_valid(st_read(raw_file, quiet = TRUE))
dist$shdist  <- as.integer(as.character(dist$REGCODE))            # boundary code = NFHS-5 district code
dist$poly_id <- seq_len(nrow(dist))

## 2. Some codes are carried by more than one polygon. Give each such polygon the
##    district to which most of the survey clusters located inside it belong.
dup_codes <- unique(dist$shdist[duplicated(dist$shdist)])
cat("Codes used by more than one polygon:", dup_codes, "\n")
a    <- readRDS(file.path(OUT_DIR, "analysis_dataset.rds"))
cl_d <- unique(as.data.table(a)[, .(DHSCLUST = hv001, shdist_cl = shdist)]); rm(a); invisible(gc())
ge <- st_read("IAGE7AFL.shp", quiet = TRUE)
ge <- ge[!(ge$LATNUM == 0 & ge$LONGNUM == 0), "DHSCLUST"]
ge <- st_transform(merge(ge, cl_d, by = "DHSCLUST"), st_crs(dist))
j    <- st_join(ge, dist[dist$shdist %in% dup_codes, "poly_id"], join = st_within)
vote <- as.data.table(st_drop_geometry(j))[!is.na(poly_id), .N, by = .(poly_id, shdist_cl)][order(poly_id, -N)]
best <- vote[, .SD[1], by = poly_id]
dist$shdist[match(best$poly_id, dist$poly_id)] <- best$shdist_cl

## 3. Checks and save
cat("Districts in the survey without a polygon:", setdiff(cl_d$shdist_cl, dist$shdist), "\n")
cat("Codes still duplicated:", unique(dist$shdist[duplicated(dist$shdist)]), "\n")
st_write(dist, "sdr_boundaries/nfhs5_districts_fixed.gpkg", delete_dsn = TRUE, quiet = TRUE)
