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
# define endpoint
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")

# we want both sentinel and landsat hls data
HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

# beaverdam reservoir bbox
bdr_box <- c(xmin = -79.827088, ymin = 37.311798, xmax = -79.811865, ymax = 37.321694)
# falling creek reservoir bbox
fcr_box <- c(xmin = -79.840037, ymin = 37.301435, xmax = -79.833651, ymax = 37.311487)
# carvins cove reservoir bbox
ccr_box <- c(xmin = -79.981728, ymin = 37.367542, xmax = -79.942552, ymax = 37.407255)
ccr_box_wgs84 <- sf::st_bbox(c(xmin = -79.981728, ymin = 37.367542, 
                               xmax = -79.942552, ymax = 37.407255), 
                             crs = "EPSG:4326")

# Convert the bounding boxes to the correct UTM projection
ccr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(ccr_box_wgs84), "EPSG:32617")
)

# set date of interest (will need to change this later to update every day)
start_date <- "2024-07-01T00:00:00Z"
end_date <- "2024-07-30T00:00:00Z"

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

# identify water using the Otsu thresholding method (https://ieeexplore.ieee.org/document/4310076)
vals <- as.vector(clean_mndwi_stars[[1]])
# scale to 0 to 1 so we can apply Otsu thresholding
vals_scaled <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
# calculate the Otsu threshold 
threshold <- otsu(matrix(vals_scaled, ncol = 1)) # use otsu to threshold
# threshold the image data to create a binary mask
binary_mask_matrix <- (clean_mndwi_stars[[1]] > threshold)
binary_stars <- st_as_stars(list(star_mask = binary_mask_matrix),
                            dimensions = st_dimensions(clean_stars))


# viz
tm_shape(shp = binary_stars) + 
  tm_raster(col = "star_mask")


test <- as(clean_stars, "Raster")
plot(test)
test_mask <- as(binary_stars, "Raster")
plot(test_mask)

test_masked <- terra::mask(test$layer, test_mask$layer)

plot(test_masked)


################################################################################
# with rast() instead of gdalcubes
################################################################################

# place items in spatial df
sf_items <- items_as_sf(items)

# retrieve granule ID for each feature
granule_id <- sapply(items$features, function(feature) feature$id)

# add as first column in sf_items
sf_items <- cbind(granule = granule_id, sf_items)

# set up function to grab bands of interest for both HLSS and HLSL
extract_asset_urls <- function(feature) {
  collection_id <- feature$collection
  if (collection_id == "HLSS30_2.0") {
    bands = c('B8A','B04','B03', 'B11', 'Fmask')
  } else if (collection_id == "HLSL30_2.0") {
    bands = c('B05','B04', 'B03', 'B11', 'Fmask')}
  sapply(bands, function(band) feature$assets[[band]]$href)
}


# retrieve asset URLs for each feature using extract_asset_urls function and transpose them to columns
asset_urls <- t(sapply(items$features, extract_asset_urls))
colnames(asset_urls) <- c('nir', 'red', 'green', 'SWIR', 'fmask') # clean up column names


sf_items <- cbind(sf_items, asset_urls)

# reset row indices
row.names(sf_items) <- NULL

## 
setGDALconfig("GDAL_HTTP_UNSAFESSL", value = "YES")
setGDALconfig("GDAL_HTTP_COOKIEFILE", value = ".rcookies")
setGDALconfig("GDAL_HTTP_COOKIEJAR", value = ".rcookies")
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", value = "EMPTY_DIR")
setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", value = "TIF")

# read an HLS scene from URL, apply the scale factor if necessary, and optionally crops and
# masks the scene based on a polygon. Requery the above GDAL configurations.
open_hls <- function(url, roi = NULL) {
  # Add VSICURL prefix
  url <- paste0('/vsicurl/', url)
  # Retrieve metadata
  meta <- describe(url)
  # Check if dataset is Quality Layer (Fmask) - no scaling this asset (int8 datatype)
  is_fmask <- any(grep("Fmask", meta))
  # Check if Scale is present in band metadata
  will_autoscale <- any(grep("Scale:", meta))
  # Read the raster
  r <- rast(url)
  # Apply Scale Factor if necessary
  if (!will_autoscale && !is_fmask){
    print(paste("No scale factor found in band metadata. Applying scale factor of 0.0001 to", basename(url)))
    r <- r * 0.0001
  }
  # Crop if roi specified
  if (!is.null(roi)){
    # project to degrees
    r_reproj <- project(r, "epsg:4326")
    # crop to roi
    e <- terra::ext(roi[1], roi[3], roi[2], roi[4])
    r <- terra::crop(r_reproj, e)
  }
  return(r)
}

swir <- open_hls(sf_items$SWIR[1], ccr_box)

ggplot() + 
  geom_spatraster(data = swir) +
  scale_fill_viridis() +
  theme_classic() 

green <- open_hls(sf_items$green[1], ccr_box)

ndwi <- (green - swir) / (green + swir)

ggplot() + 
  geom_spatraster(data = ndwi) +
  scale_fill_viridis() +
  theme_classic() 
# identify water using the Otsu thresholding method
vals <- values(ndwi, mat = F) # get ndwi values
vals <- vals[!is.na(vals)] # remove any nas
vals_scaled <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
threshold <- otsu(matrix(vals_scaled, ncol = 1)) # use otsu to threshold
water_mask <- ndwi > threshold # identify values over threshold
#classify and plot
water_mask01 <- classify(water_mask, rcl = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))
plot(water_mask01, main = "Water Mask (Otsu Threshold)")
plot(water_mask)
ggplot() + 
  geom_spatraster(data = water_mask) +
  scale_fill_viridis() +
  theme_classic() 

# extract value at location
location <- data.frame(x = 37.315348, y = -79.819365)
pt <- data.frame(lon = -79.819365, lat = 37.315348)
nir_val <- terra::extract(nir_reproj, pt)




