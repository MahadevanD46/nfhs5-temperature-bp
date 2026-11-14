library(ecmwfr)
wf_set_key(key = "PASTE_YOUR_CDS_API_TOKEN")

req_norm <- list(
  dataset_short_name = "reanalysis-era5-land-monthly-means",
  product_type    = "monthly_averaged_reanalysis",
  variable        = "2m_temperature",
  year            = as.character(1991:2018),
  month           = sprintf("%02d", 1:12),
  time            = "00:00",
  data_format     = "netcdf",
  download_format = "unarchived",
  area            = c(38, 68, 6, 98),
  target          = "t2m_monthly_1991_2018.nc")

wf_request(req_norm, path = "climate", transfer = TRUE)

library(terra)
rm_ <- rast("climate/t2m_monthly_1991_2018.nc")
nlyr(rm_)            # should be 336 (28 years x 12 months)
range(time(rm_))     # 1991-01-01 to 2018-12-01

library(terra); library(data.table)
rm_ <- rast("climate/t2m_monthly_1991_2018.nc")
mon <- as.integer(format(time(rm_), "%m"))
vm  <- as.matrix(terra::extract(rm_, pts, ID = FALSE)) - 273.15     # clusters x 336 months
clim12 <- sapply(1:12, function(m) rowMeans(vm[, mon == m, drop = FALSE], na.rm = TRUE))  # clusters x 12

# Harmonic regression (2 harmonics) through mid-month points -> daily normal for doy 1..366
H <- function(doy) cbind(1, sin(2*pi*doy/365.25), cos(2*pi*doy/365.25),
                         sin(4*pi*doy/365.25), cos(4*pi*doy/365.25))
mid <- c(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
B   <- solve(crossprod(H(mid)), t(H(mid)) %*% t(clim12))           # 5 x clusters
norm_daily <- H(1:366) %*% B                                        # 366 x clusters

annual_mean <- data.table(DHSCLUST = ge$DHSCLUST, t_annual = rowMeans(clim12))
saveRDS(list(norm_daily = norm_daily, annual_mean = annual_mean, ids = ge$DHSCLUST),
        "output/normals.rds")
