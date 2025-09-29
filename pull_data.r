# Code started by Molly Stroud on 9/29/25

# load in packages
require(pacman)
p_load('rmarkdown','earthdatalogin', 'rstac','imager','lubridate','xts',
       'dygraphs','leaflet','terra')

# define endpoint
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")

HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

# set bounding box of beaverdam reservoir
bdr_box <- c(xmin = -79.827088, ymin = 37.311798, xmax = -79.811865, ymax = 37.321694)

# set date of interest
roi_datetime <- '2021-08-01T00:00:00Z/2021-08-15T23:59:59Z' 

# see what's available
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = bdr_box,
              datetime = roi_datetime,
              limit = 100) %>%
  post_request()
print(items)

# select an image and plot it
browse_image_url <- items$features[[3]]$assets$browse$href
browse <- imager::load.image(browse_image_url)
#plot(browse)

# place items in spatial df
sf_items <- items_as_sf(items)

# Retrieve Granule ID for each feature
granule_id <- sapply(items$features, function(feature) feature$id)
# Add as first column in sf_items
sf_items <- cbind(granule = granule_id, sf_items)
View(sf_items)

# select bands, filter clouds
# Define a function to extract asset urls for selected bands
# # This also includes a check to ensure the correct bands are extracted
# # # depending on the collection (HLSL30 or HLSS30)
extract_asset_urls <- function(feature) {
  collection_id <- feature$collection
  if (collection_id == "HLSS30_2.0") {
    bands = c('B8A','B04','Fmask')
  } else if (collection_id == "HLSL30_2.0") {
    bands = c('B05','B04','Fmask')}
  sapply(bands, function(band) feature$assets[[band]]$href)
}
# Retrieve Asset URLs for each feature using our extract_asset_urls function and transpose them to columns
asset_urls <- t(sapply(items$features, extract_asset_urls))
View(asset_urls)
colnames(asset_urls) <- c('nir', 'red', 'fmask')
sf_items <- cbind(sf_items, asset_urls)
View(sf_items)

# Filter based on cloud cover
sf_items <- sf_items[sf_items$eo.cloud_cover < 30,]
# Reset Row Indices
row.names(sf_items) <- NULL
View(sf_items)


setGDALconfig("GDAL_HTTP_UNSAFESSL", value = "YES")
setGDALconfig("GDAL_HTTP_COOKIEFILE", value = ".rcookies")
setGDALconfig("GDAL_HTTP_COOKIEJAR", value = ".rcookies")
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", value = "EMPTY_DIR")
setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", value = "TIF")

# This function reads an HLS scene from a URL, applies the scale factor if necessary, and optionally crops and
# masks the scene based on a polygon. It requries the above GDAL configurations and a .netrc file. A .netrc
# can be created by running `earthdatalogin::edl_netrc()`.
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
red <- open_hls(sf_items$red[3])
plot(red)

nir <- open_hls(sf_items$nir[3])
plot(nir)

# make a stacked raster
stacked <- c(red, nir)
stacked_reproj <- project(stacked,"epsg:4326")
stacked_crop <- crop(stacked_reproj, bdr_box)

plot(stacked_crop)

# make stacks
red_stack <- lapply(sf_items$red, open_hls)
nir_stack <- lapply(sf_items$nir, open_hls)
fmask_stack <- lapply(sf_items$fmask, open_hls)

################################################################################
# build mask to grab water
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


