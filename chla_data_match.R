################################################################################
# Code started by Molly Stroud on 10/10/25
################################################################################

# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf',
       'fuzzyjoin')
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

# convert the bounding boxes to the correct UTM projection
bvr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.827088, 
                                               ymin = 37.311798, 
                                               xmax = -79.811865, 
                                               ymax = 37.321694), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)
fcr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.840037, 
                                               ymin = 37.301435, 
                                               xmax = -79.833651, 
                                               ymax = 37.311487), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)
ccr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.981728, 
                                               ymin = 37.367522, 
                                               xmax = -79.942552, 
                                               ymax = 37.407255), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)

# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")


# function to identify all imagery within start and end dates
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
    ext_query("eo:cloud_cover" < 20) %>% # filter for cloud cover
    post_request() %>%
    items_fetch()
  dates <- sapply(items$features, function(x) x$properties$datetime)
}

# call function to get dates
ccr_hls_datetime <- get_dates(ccr_box, ccr_dates)
fcr_hls_datetime <- get_dates(fcr_box, fcr_dates)
bvr_hls_datetime <- get_dates(bvr_box, bvr_dates)

# now fuzzy match dates with in situ dates
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
    max_dist = 1, # day difference
    distance_col = "day_diff"
  )
}

ccr_joined <- fuzzy_join_dates(ccr_hls_datetime, ccr_dates)
fcr_joined <- fuzzy_join_dates(fcr_hls_datetime, fcr_dates)
bvr_joined <- fuzzy_join_dates(bvr_hls_datetime, bvr_dates)

# now join back with full chla data. possibly make this a function or add to above function
# ccr
chla_ccr$DateTime <- as.Date(chla_ccr$DateTime)
ccr_chla_joined <- chla_ccr[chla_ccr$DateTime %in% ccr_joined$chla_dates,]
ccr_chla_joined <- left_join(ccr_chla_joined, ccr_joined, 
                             by = c("DateTime" = "chla_dates"))
# fcr
chla_fcr$DateTime <- as.Date(chla_fcr$DateTime)
fcr_chla_joined <- chla_fcr[chla_fcr$DateTime %in% fcr_joined$chla_dates,]
fcr_chla_joined <- left_join(fcr_chla_joined, fcr_joined, 
                             by = c("DateTime" = "chla_dates"))
#bvr
chla_bvr$DateTime <- as.Date(chla_bvr$DateTime)
bvr_chla_joined <- chla_bvr[chla_bvr$DateTime %in% bvr_joined$chla_dates,]
bvr_chla_joined <- left_join(bvr_chla_joined, bvr_joined, 
                             by = c("DateTime" = "chla_dates"))


################################################################################
# next: download data and extract points of interest using gdalcubes extract_geom
################################################################################
# ok maybe 
start_date <- paste0(bvr_chla_joined$HLS_dates[1], "T00:00:00Z")
end_date <- paste0(bvr_chla_joined$HLS_dates[30], "T00:00:00Z")
# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = bvr_box,
              datetime = paste(start_date, end_date, sep="/"),
              limit = 500) %>%
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

# sentinel collection
col_S30 <- stac_image_collection(
  ss_items,
  asset_names = c("B02", "B03", "B04", "B8A"),
  url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
)
# landsat collection
col_L30 <- stac_image_collection(
  sl_items,
  asset_names = c("B02", "B03", "B04", "B05"),
  url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
)

# define the cube space
cube <- cube_view(srs ="EPSG:32617",
                  extent = list(t0 = start_date, 
                                t1 = end_date,
                                left = bvr_box_utm[1], 
                                right = bvr_box_utm[3],
                                top = bvr_box_utm[4], 
                                bottom = bvr_box_utm[2]),
                  dx = 30, # 30 m resolution
                  dy = 30, 
                  dt = "P1D",
                  aggregation = "median", 
                  resampling = "average")

# HLSS
data_S30 <- raster_cube(image_collection = col_S30, 
                        view = cube)
data_S30 <- rename_bands(data_S30, B02 = "blue", B03 = "green", B04 = "red",
                         B8A = "NIR")
#HLSL
data_L30 <- raster_cube(image_collection = col_L30, 
                        view = cube)
data_L30 <- rename_bands(data_L30, B02 = "blue", B03 = "green", B04 = "red",
                         B05 = "NIR")

# get unique lat/longs
points <- data.frame(unique(cbind(bvr_chla_joined$Longitude, bvr_chla_joined$Latitude)))
points <- mutate(points, FID = row_number())
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- st_as_sf(x = points,                         
         coords = c("X1", "X2"),
         crs = projcrs)
test <- extract_geom(data_S30, points, drop_geom = F)
test


with_coords <- left_join(test, points, by = "FID")

