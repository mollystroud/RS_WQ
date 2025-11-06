################################################################################
# Code started by Molly Stroud on 9/29/25
################################################################################
# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf')
# likely unused: 'imager', 'lubridate','xts','dygraphs', 'leaflet', 'cubelyr'
################################################################################
## the below code is designed to pull landsat thermal imagery over a specified 
# area and estimate temperature over the reservoir
################################################################################
# get bboxes
source("roanokeres_bbox.R")

# define stac url
ls = stac("https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2")
ls = stac ("https://planetarycomputer.microsoft.com/api/stac/v1")
# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")


items <- ls %>%
  stac_search(collections = "landsat-c2-l2",
              bbox = fcr_box,
              datetime = paste("2020-01-01T00:00:00Z", "2020-03-01T00:00:00Z", sep="/"),
              limit = 500) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request() %>%
  items_sign(sign_fn = sign_planetary_computer()) %>%
  items_fetch()
items
  
items$features[[3]]$properties$`proj:epsg`
items$features <- Filter(
  function(f) f$properties$platform %in% c("landsat-8", "landsat-9"),
  items$features
)

cube <- cube_view(srs ="EPSG:32617",
                  extent = list(t0 = "2019-01-01T00:00:00Z", 
                                t1 = "2019-03-01T00:00:00Z",
                                left = fcr_box_utm[1], 
                                right = fcr_box_utm[3],
                                top = fcr_box_utm[4], 
                                bottom = fcr_box_utm[2]),
                  dx = 30, # 30 m resolution
                  dy = 30, 
                  dt = "P1D",
                  aggregation = "median", 
                  resampling = "average")

col <- stac_image_collection(items$features,
  asset_names = c("red", "green", "blue", "lwir11"),
  url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
)

# make raster cube
data <- raster_cube(image_collection = col, 
                    view = cube)
data <- rename_bands(data, red = "red", green = "green", blue = "blue",
                     lwir11 = "thermal")

test <- st_as_stars(data)
test
plot(test)



################################################################################
# Function to create data cube object with set dates and bbox
################################################################################
get_lst <- function(bbox, bbox_utm, start_date, end_date) {
  # grab items within dates of interest
  items <- ls %>%
    stac_search(collections = "landsat-c2-l2",
                bbox = bbox,
                datetime = paste(start_date, end_date, sep="/"),
                limit = 500) %>%
    ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
    post_request() %>%
    items_sign(sign_fn = sign_planetary_computer()) %>%
    items_fetch()
  # filter out non LS8/9
  items$features <- Filter(
    function(f) f$properties$platform %in% c("landsat-8", "landsat-9"),
    items$features
  )

  # define the cube space
  cube <- cube_view(srs ="EPSG:32617",
                    extent = list(t0 = start_date, 
                                  t1 = end_date,
                                  left = bbox_utm[1], 
                                  right = bbox_utm[3],
                                  top = bbox_utm[4], 
                                  bottom = bbox_utm[2]),
                    dx = 30, # 30 m resolution
                    dy = 30, 
                    dt = "P1D",
                    aggregation = "median", 
                    resampling = "average")
    # create stac image collection
    col <- stac_image_collection(items$features,
                               asset_names = c("red", "green", "blue", "lwir11"),
                               url_fun = function(url) paste0("/vsicurl/", url))  # helps GDAL access
    # make raster cube
    data <- raster_cube(image_collection = col, 
                        view = cube)
    data <- rename_bands(data, red = "red", green = "green", blue = "blue",
                         lwir11 = "thermal")
    # make stars obj
    ls_stars <- st_as_stars(data)
    # remove empty dates
    arr <- ls_stars[[1]] # extract raw array (x, y, time)
    non_na_counts <- apply(arr, 3, function(slice) sum(!is.na(slice))) # count non-NA pixels for each time
    valid_idx <- which(non_na_counts > 0) # indices of slices that have at least one real value
    # build cleaned object by stacking only valid slices
    slices <- lapply(valid_idx, function(i) ls_stars[,,, i, drop = FALSE])
    clean_ls_stars <- do.call(c, c(slices, along = "time"))
  return(clean_ls_stars)
}

# set date of interest
start_date <- "2022-01-01T00:00:00Z"
end_date <- "2022-03-30T00:00:00Z"
data_lst <- get_lst(fcr_box, fcr_box_utm, start_date, end_date)

tm_shape(shp = data_lst) + 
  tm_rgb(r = "red", g = "green", b = "blue")

thermal <- data_lst["thermal"]
# convert to kelvin
thermal_celsius <- (thermal * 0.00341802) + 149.0 - 273.15
names(thermal_celsius) <- "temp_C"
plot(thermal_celsius, main = "Land Surface Temperature (°C)")
tm_shape(shp = thermal_celsius) +
  tm_raster(style = "cont")

