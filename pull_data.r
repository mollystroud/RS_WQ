################################################################################
# Code started by Molly Stroud on 9/29/25
################################################################################
# load in packages
require(pacman)
p_load('rmarkdown','earthdatalogin', 'rstac','imager','lubridate','xts',
       'dygraphs','leaflet','terra', 'stars', 'ggplot2', 'tidyterra', 'viridis',
       'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'cubelyr')
################################################################################
## the below code is designed to pull HLS imagery over a specified area and 
# mask out any land
################################################################################
# define stac url
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")

# both sentinel and landsat hls data
HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

# beaverdam reservoir bbox
bdr_box <- c(xmin = -79.827088, 
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
             ymin = 37.367542, 
             xmax = -79.942552, 
             ymax = 37.407255)

# convert the bounding boxes to the correct UTM projection
bdr_box_utm <- sf::st_bbox(
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
                                               ymin = 37.367542, 
                                               xmax = -79.942552, 
                                               ymax = 37.407255), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)

# set date of interest (will need to change this later to update every day)
start_date <- "2025-09-01T00:00:00Z"
end_date <- "2025-09-30T00:00:00Z"

# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = ccr_box,
              datetime = paste(start_date, end_date, sep="/"),
              limit = 100) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request()
items

# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")

# manually add projection info into each STAC item before passing to gdalcubes
for (i in seq_along(items$features)) {
  items$features[[i]]$properties$`proj:epsg` <- 32617  # force correct projection
}

# pass to gdalcubes
# need to separate HLSL and HLSS data since they have different band names
ss_items <- items$features[sapply(items$features, function(f) f$collection == "HLSS30_2.0")]
sl_items <- items$features[sapply(items$features, function(f) f$collection == "HLSL30_2.0")]
# sentinel collection
col_S30 <- stac_image_collection(
  ss_items,
  asset_names = c("B02", "B03", "B04", "B8A", "B11"),
  url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
)
# landsat collection
col_L30 <- stac_image_collection(
  sl_items,
  asset_names = c("B02", "B03", "B04", "B05", "B06"),
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

# make raster cube
data_S30 <- raster_cube(image_collection = col_S30, 
                     view = cube)
data_S30 <- rename_bands(data_S30, B02 = "blue", B03 = "green", B04 = "red",
                         B8A = "NIR", B11 = "SWIR")

data_L30 <- raster_cube(image_collection = col_L30, 
                        view = cube)
data_L30 <- rename_bands(data_L30, B02 = "blue", B03 = "green", B04 = "red",
                         B05 = "NIR", B06 = "SWIR")

# test plot
#plot(data, rgb=3:1, zlim=c(0,1200))

# calculate mNDWI
mndwi_S30 <- data_S30 |>
  select_bands(bands = c("green", "SWIR")) |>
  apply_pixel(expr = "(green-SWIR)/(green+SWIR)", names = "mNDWI")

mndwi_L30 <- data_L30 |>
  select_bands(bands = c("green", "SWIR")) |>
  apply_pixel(expr = "(green-SWIR)/(green+SWIR)", names = "mNDWI")

# open as stars objects and merge
mndwi_stars_S30 <- st_as_stars(mndwi_S30)
mndwi_stars_L30 <- st_as_stars(mndwi_L30)
mndwi_stars <- c(mndwi_stars_S30, mndwi_stars_L30, along = 3)

# remove dates with no imagery
# extract raw array (x, y, time)
arr <- mndwi_stars[[1]]

# count non-NA pixels for each time slice
non_na_counts <- apply(arr, 3, function(slice) sum(!is.na(slice)))

# indices of slices that have at least one real value
valid_idx <- which(non_na_counts > 0)

# build cleaned object by stacking only valid slices
slices <- lapply(valid_idx, function(i) mndwi_stars[,,, i, drop = FALSE])
clean_mndwi_stars <- do.call(c, c(slices, along = "time"))

# viz
tm_shape(shp = clean_mndwi_stars) + 
  tm_raster(col = "mNDWI")

# identify water using the Otsu thresholding method for each image 
# (https://ieeexplore.ieee.org/document/4310076)
n_dates <- dim(clean_mndwi_stars)["time"] # get dates
time_labels <- st_get_dimension_values(clean_mndwi_stars, "time")

# get spatial dims
nx <- dim(clean_mndwi_stars)["x"]
ny <- dim(clean_mndwi_stars)["y"]

# set up for loop
mask_array <- array(NA, dim = c(nx, ny, n_dates))
thresholds <- numeric(n_dates)
# for each date, scale mndwi values (0-1) and calculate Otsu threshold
for (i in 1:n_dates) {
  slice <- clean_mndwi_stars[,,,i, drop = TRUE]
  vals <- as.vector(slice[[1]])
  
  vals_scaled <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
  threshold_i <- otsu(matrix(vals_scaled, ncol = 1))
  thresholds[i] <- threshold_i
  
  mask_array[,,i] <- slice[[1]] > threshold_i
}

# convert array back to stars
binary_stars <- st_as_stars(mask_array)

# attach dimensions (x, y from original, time from labels)
binary_stars <- st_set_dimensions(binary_stars, "x", 
                                  values = st_get_dimension_values(clean_mndwi_stars, "x"))
binary_stars <- st_set_dimensions(binary_stars, "y", 
                                  values = st_get_dimension_values(clean_mndwi_stars, "y"))
binary_stars <- st_set_dimensions(binary_stars, "time", 
                                  values = time_labels)

names(binary_stars) <- "watermask"

# viz
tm_shape(shp = binary_stars) + 
  tm_raster(col = "watermask")





test <- as(clean_stars, "Raster")
plot(test)
test_mask <- as(binary_stars, "Raster")
plot(test_mask)

test_masked <- terra::mask(test$layer, test_mask$layer)

plot(test_masked)

