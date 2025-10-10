################################################################################
# Code started by Molly Stroud on 10/10/25
################################################################################

# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf')
################################################################################
# the below code is designed to match up in situ chl-a data over CCR, FCR, and
# BVR with HLS data
################################################################################
# read in filtered chl-a data from reservoirs
chla <- read_csv("filt-chla_2014_2024.csv")
# only surface measurements
chla <- chla[chla$Depth_m <= 0.1,]
# match with lat/long
sitelocs <- read_csv("site_descriptions.csv")
chla <- left_join(chla, sitelocs[2:5], by = "Site")
# only ccr, bvr, fcr
chla_ccr <- chla[chla$Reservoir == "CCR",]
chla_fcr <- chla[chla$Reservoir == "FCR",]
chla_bvr <- chla[chla$Reservoir == "BVR",]
# get unique dates
ccr_dates <- unique(as.Date(chla_ccr$DateTime))
fcr_dates <- unique(as.Date(chla_fcr$DateTime))
bvr_dates <- unique(as.Date(chla_bvr$DateTime))

# define stac url
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")

# both sentinel and landsat hls data
HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

# beaverdam reservoir bbox
bvr_box <- c(xmin = -79.827088, 
             ymin = 37.311798, 
             xmax = -79.811865, 
             ymax = 37.321694)
# falling creek reservoir bbox
fcr_box <- c(xmin = -79.840037, 
             ymin = 37.301435, 
             xmax = -79.833651, 
             ymax = 37.311487)
# carvins cove reservoir bbox
ccr_box <- c(xmin = -79.981728, 
             ymin = 37.367522, 
             xmax = -79.942552, 
             ymax = 37.407255)

# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")


# grab items within dates of interest
get_dates <- function(bbox, dates) {
  # get start and end dates
  start_date <- paste0(dates[1], "T00:00:00Z")
  end_date <- paste0(dates[length(dates)], "T00:00:00Z")
  # search
  items <- s %>%
    stac_search(collections = HLS_col,
                bbox = bbox,
                datetime = paste(start_date, end_date, sep="/"),
                limit = 500) %>%
    ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
    post_request() %>%
    items_fetch()
  dates <- sapply(items$features, function(x) x$properties$datetime)
}

# call function to get dates
ccr_hls_datetime <- get_dates(ccr_box, ccr_dates)
fcr_hls_datetime <- get_dates(fcr_box, fcr_dates)
bvr_hls_datetime <- get_dates(bvr_box, bvr_dates)

# now fuzzy match dates with in situ dates
# first, lose the time
# ccr
ccr_hls_dates <- t(data.frame(lapply(ccr_hls_datetime, as.Date)))
row.names(ccr_hls_dates) <- NULL
colnames(ccr_hls_dates) <- c("Date")
# fcr
fcr_hls_dates <- t(data.frame(lapply(fcr_hls_datetime, as.Date)))
row.names(fcr_hls_dates) <- NULL
colnames(fcr_hls_dates) <- c("Date")
# bvr
bvr_hls_dates <- t(data.frame(lapply(bvr_hls_datetime, as.Date)))
row.names(bvr_hls_dates) <- NULL
colnames(bvr_hls_dates) <- c("Date")





