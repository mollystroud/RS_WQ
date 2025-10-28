################################################################################
# Code started by Molly Stroud on 10/21/25
################################################################################
# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf',
       'fuzzyjoin', 'patchwork')

################################################################################
# the below code will estimate chl-a values over the Upper Yarra Reservoir
################################################################################
################################################################################
# read in WQ from Upper Yarra
################################################################################
yarra_wq <- read_csv("YARR-targets-insitu.csv")
yarra_chla <- yarra_wq[yarra_wq$variable == "PHY_TCHLA",]
yarra_chla <- yarra_wq[yarra_wq$variable == "turbidity",]

yarra_chla <- yarra_chla[yarra_chla$depth <= 0.1,]
yarra_dates <- unique(as.Date(yarra_chla$datetime[yarra_chla$datetime > "2013-04-10"]))


################################################################################
# set up to pull from stac
################################################################################
# define stac url
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")
# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")

# both sentinel and landsat hls data
HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

# upper yarra reservoir bbox
yarra_box <- c(xmin = 145.886, 
             ymin = -37.723, 
             xmax = 145.983, 
             ymax = -37.642)
# convert the bounding boxes to the correct UTM projection
yarra_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = 145.886, 
                                               ymin = -37.723, 
                                               xmax = 145.983, 
                                               ymax = -37.642), 
                                             crs = "EPSG:4326")), "EPSG:32617"))
# set date of interest
start_date <- "2013-04-12T00:00:00Z"
end_date <- "2025-06-18T00:00:00Z"
# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = yarra_box,
              datetime = paste(start_date, end_date, sep="/"),
              limit = 100) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request()
items

# manually add projection info into each STAC item before passing to gdalcubes
# for some reason HLS data doesn't have this automatically included?
for (i in seq_along(items$features)) {
  items$features[[i]]$properties$`proj:epsg` <- 32617  # force correct projection
}

# pass to gdalcubes
# need to separate HLSL and HLSS data since they have different band configs
ss_items <- items$features[sapply(items$features, function(f) f$collection == "HLSS30_2.0")]
sl_items <- items$features[sapply(items$features, function(f) f$collection == "HLSL30_2.0")]

################################################################################
# function to identify all imagery within start and end dates
################################################################################
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
    ext_query("eo:cloud_cover" < 10) %>% # filter for cloud cover
    post_request() %>%
    items_fetch()
  dates <- sapply(items$features, function(x) x$properties$datetime)
}

# call function to get dates
yarra_hls_datetime <- get_dates(yarra_box, (yarra_dates)) # need to order in ascending chrono, rev if needed



################################################################################
# function to fuzzy match dates with in situ dates
################################################################################
fuzzy_join_dates <- function(hls_datetime, chla_dates){
  # put into dataframe
  hls_dates <- data.frame(t(data.frame(lapply(hls_datetime, as.Date))))
  row.names(hls_dates) <- NULL # remove row names
  colnames(hls_dates) <- c("HLS_dates") # name column
  hls_dates$HLS_dates <- as.Date(hls_dates$HLS_dates) # make Date
  # join together
  chla_dates <- data.frame(chla_dates)
  joined <- difference_inner_join(
    hls_dates, chla_dates,
    by = c("HLS_dates" = colnames(chla_dates)),
    max_dist = 2, # day difference
    distance_col = "day_diff"
  )
}

yarra_joined <- fuzzy_join_dates(yarra_hls_datetime, yarra_dates)

# now join back with full chla data
yarra_chla$datetime <- as.Date(yarra_chla$datetime)
yarra_chla_joined <- yarra_chla[yarra_chla$datetime %in% yarra_joined$chla_dates,]
yarra_chla_joined <- left_join(yarra_chla_joined, yarra_joined, 
                             by = c("datetime" = "chla_dates"))




