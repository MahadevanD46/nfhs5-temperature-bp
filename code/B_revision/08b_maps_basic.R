# 08b_maps_basic.R — plain maps of the state and district results (default styling)
# State outline: Survey of India boundaries from the imdR package.
source("code/B_revision/00_config.R")
suppressPackageStartupMessages({ library(data.table); library(sf); library(ggplot2) })
sf_use_s2(FALSE)
st <- st_transform(st_make_valid(st_as_sf(imdR::india_states)), 7755)
name_col <- names(st)[sapply(st_drop_geometry(st), function(x) is.character(x) && any(grepl("punjab", x, ignore.case = TRUE)))][1]
key <- function(x) {
  x <- gsub("[^a-z]", "", gsub("^the |nct of |national capital territory of ", "", gsub("&", "and", tolower(x))))
  x[grepl("andaman", x)] <- "andaman"; x[grepl("dadra|daman", x)] <- "dnhdd"
  x[grepl("jammu", x)] <- "jammuandkashmir"; x[grepl("delhi", x)] <- "delhi"; x
}
st$key <- key(st[[name_col]])

## state maps: temperature difference and change in prevalence
t3 <- fread(file.path(OUT_DIR, "table3_bootstrap.csv"))
a  <- load_a(cols = c("hv024", "t0", "t_annual", "wt"))
sh <- a[!is.na(t_annual) & !is.na(t0), .(w_shift = weighted.mean(t0 - t_annual, wt)), by = hv024]; rm(a)
tab <- merge(t3[, .(hv024, state, change)], sh, by = "hv024")[hv024 != 37]      # Ladakh not shown
tab[, key := key(state)]
mp <- merge(st, tab, by = "key", all.x = TRUE)
for (v in c("w_shift", "change"))
  ggsave(file.path(FIG_DIR, paste0("map_state_", v, ".png")),
         ggplot(mp) + geom_sf(aes(fill = .data[[v]])) + scale_fill_gradient2(), width = 6, height = 6, dpi = 150)

## district LISA classes
dp <- fread(file.path(OUT_DIR, "lisa_point.csv"))
d  <- st_read(DIST_SHP, quiet = TRUE); d$shdist <- as.integer(as.character(d[[DIST_ID_FIELD]]))
md <- merge(st_transform(d[, "shdist"], 7755), dp[, .(shdist, class_obs, class_std)], by = "shdist", all.x = TRUE)
for (v in c("class_obs", "class_std"))
  ggsave(file.path(FIG_DIR, paste0("map_lisa_", v, ".png")),
         ggplot() + geom_sf(data = md, aes(fill = .data[[v]]), linewidth = 0.05) + geom_sf(data = st, fill = NA),
         width = 6, height = 6, dpi = 150)
log_msg("08b done")
