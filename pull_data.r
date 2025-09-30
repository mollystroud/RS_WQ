# Code started by Molly Stroud on 9/29/25

# load in packages
require(pacman)
p_load('rmarkdown','earthdatalogin', 'rstac','imager','lubridate','xts',
       'dygraphs','leaflet','terra', 'gdalcubes', 'stars', 'ggplot2', 'tidyterra')

################################################################################
#### following nasa tutorial
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
ccr_box <- c(xmin = -79.981728, ymin = 37.367572, xmax = -79.942552, ymax = 37.407255)

# set date of interest (will need to change this later to update every day)
datetime <- '2021-08-01T00:00:00Z/2021-08-10T23:59:59Z' 

# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = bdr_box,
              datetime = datetime,
              limit = 100) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request()
items


# set up function to grab bands of interest for both HLSS and HLSL
extract_asset_urls <- function(feature) {
  collection_id <- feature$collection
  if (collection_id == "HLSS30_2.0") {
    bands = c('B8A','B04','Fmask')
  } else if (collection_id == "HLSL30_2.0") {
    bands = c('B05','B04','Fmask')}
  sapply(bands, function(band) feature$assets[[band]]$href)
}

# place items in spatial df
sf_items <- items_as_sf(items)

# retrieve granule ID for each feature
granule_id <- sapply(items$features, function(feature) feature$id)

# add as first column in sf_items
sf_items <- cbind(granule = granule_id, sf_items)

# retrieve asset URLs for each feature using extract_asset_urls function and transpose them to columns
asset_urls <- t(sapply(items$features, extract_asset_urls))
colnames(asset_urls) <- c('nir', 'red', 'fmask') # clean up column names
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
    # Reproject roi to match crs of r
    roi_reproj <- project(roi, crs(r))
    r <- mask(crop(r, roi_reproj), roi_reproj)
  }
  return(r)
}

# Test opening and crop
nir <- open_hls(sf_items$nir[2])
nir_reproj <- project(nir, "epsg:4326")
e <- terra::ext(-79.827088, -79.811865, 37.31, 37.324)
nir_crop <- terra::crop(nir_reproj, e)

ggplot() + 
  geom_spatraster(data = nir_crop) +
  theme_classic() 

# extract value at location
location <- data.frame(x = 37.315348, y = -79.819365)
pt <- data.frame(lon = -79.819365, lat = 37.315348)
nir_val <- terra::extract(nir_reproj, pt)




################################################################################
# build mask to exclude non-water pixels
################################################################################
build_mask <- function(fmask, selected_bit_nums){
  # Create a mask of all zeros
  mask <- rast(fmask, vals=0)
  for (b in selected_bit_nums){
    # Apply Bitwise AND to fmask values and selected bit numbers
    mask_temp <- app(fmask, function(x) bitwAnd(x, bitwShiftL(1,b)) >0)
    # Update Mask to maintain only 1 layer with bitwise OR
    mask <- mask | mask_temp
  }
  return(mask)
}
selected_bit_nums <- c(1,2,3,4)

qmask <- build_mask(fmask, selected_bit_nums)
plot(qmask)

# Create List of Masks
qmask_stack <- lapply(fmask_stack, build_mask, selected_bit_nums=selected_bit_nums)

# Apply Mask using NA Values
masked <- mapply(function(x, y) {
  mask(x, y, maskvalue = TRUE, updatevalue = NA)
}, red_stack, qmask_stack, SIMPLIFY = FALSE)


masked <- rast(masked)

plot(masked[[2]])




