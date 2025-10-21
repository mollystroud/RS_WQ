################################################################################
# Code started by Molly Stroud on 10/21/25
################################################################################
# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf',
       'fuzzyjoin', 'patchwork')

################################################################################
# the below code will estimate chl-a values using the algorithm created in
# chla_algorithm.R, either at a point location or over the entire reservoir,
# over the specified dates
################################################################################
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

# set date of interest (will need to change this to update every day)
start_date <- "2024-01-01T00:00:00Z"
end_date <- "2024-12-31T00:00:00Z"

# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")

# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = ccr_box,
              datetime = paste(start_date, end_date, sep="/"),
              limit = 100) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request()
items

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
                                left = ccr_box_utm[1], 
                                right = ccr_box_utm[3],
                                top = ccr_box_utm[4], 
                                bottom = ccr_box_utm[2]),
                  dx = 30, # 30 m resolution
                  dy = 30, 
                  dt = "P1D",
                  aggregation = "median", 
                  resampling = "average")


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
points <- data.frame(unique(cbind(-79.9582, 37.3700)))
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # project
# make sf for extract_geom
points_crs <- st_as_sf(x = points,                         
                       coords = c("X1", "X2"),
                       crs = projcrs)
# extract band values at geometry points
# HLSS
band_vals_S30 <- extract_geom(data_S30, points_crs)
# HLSL
band_vals_L30 <- extract_geom(data_L30, points_crs)
# join
band_vals_HLS <- data.frame(rbind(band_vals_S30, band_vals_L30))
band_vals_HLS$time <- as.Date(band_vals_HLS$time) # make Date

# now estimate chl-a from GAM model in chla_algorithm.R
band_vals_HLS <- band_vals_HLS[band_vals_HLS$blue > 0 & band_vals_HLS$green > 0 &
                                 band_vals_HLS$red > 0 & band_vals_HLS$NIR > 0,]
band_vals_HLS$prediction <- exp(predict(gam_fit, newdata=band_vals_HLS)) - 0.01

# plot alongside EXO
chla_ccr_EXO <- read_csv("ccre-waterquality_2021_2024.csv")
chla_ccr_EXO <- chla_ccr_EXO[!is.na(chla_ccr_EXO$EXOChla_ugL_1),] %>%
  select(DateTime, EXOChla_ugL_1)
chla_ccr_EXO$DateTime <- as.Date(chla_ccr_EXO$DateTime)
chla_ccr_EXO <- chla_ccr_EXO %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1))

ggplot() +
  geom_line(data = band_vals_HLS, aes(x = time, y = prediction, color = 'HLS')) +
  geom_line(data = chla_ccr_EXO[chla_ccr_EXO$DateTime >= "2024-01-01",],
            aes(x = DateTime, y = mean_chla, color = 'EXO')) +
  theme_classic() +
  scale_color_manual(
    values = c("HLS" = "#065A82", "EXO" = "#FF7700")) +
  labs(x = element_blank(), y = 'Chla (uGL)', color = element_blank())

# now look for some spatial variation



