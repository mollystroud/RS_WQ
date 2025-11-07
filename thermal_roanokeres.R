################################################################################
# Code started by Molly Stroud on 9/29/25
################################################################################
# load in packages
require(pacman)
p_load('rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 'viridis', 
       'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf')
################################################################################
## the below code is designed to pull landsat thermal imagery over a specified 
# area and estimate temperature over the reservoir
################################################################################
# get bboxes
source("roanokeres_bbox.R")
# define stac url
ls = stac ("https://planetarycomputer.microsoft.com/api/stac/v1")

################################################################################
# Function to create thermal stars object with specified dates and bbox
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
                               asset_names = c("lwir11", "qa_pixel"),
                               url_fun = function(url) paste0("/vsicurl/", url))  # helps GDAL access
    # make raster cube
    data <- raster_cube(image_collection = col, 
                        view = cube)
    data <- rename_bands(data, lwir11 = "thermal", qa_pixel = "QA")
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

################################################################################
# function to get only water pixels
################################################################################
water_mask <- function(thermal_data) {
  # make array
  thermal_arr <- thermal_data$thermal
  qa_arr <- thermal_data$QA
  # get only water with bitwise flags
  water <- (bitwAnd(qa_arr, bitwShiftL(1, 7)) == 0)
  # get water or NA pixels
  thermal_arr[water | is.na(qa_arr)] <- NA_real_
  # convert to celsius
  thermal_C <- (thermal_arr * 0.00341802 + 149) - 273.15
  # make stars object again
  thermal_masked_vec <- data_lst["thermal"]   # copy dimensions
  thermal_masked_vec$thermal_C <- thermal_C
  return(thermal_masked_vec)
}

################################################################################
# next, we will extract just points in middle of water body
# hopefully this will avoid any mixed pixel issues due to the thermal band
# resampling to 30 m
################################################################################
fcr_points <- data.frame(x = c(603010, 603050), y = c(4129520, 4129150))
fcr_points <- st_as_sf(
  fcr_point,
  coords = c("x", "y"),   # change to your column names
  crs = st_crs(thermal_masked_vec)    # match raster CRS!
)
# bvr
bvr_points <- data.frame(x = c(604600), y = c(4130460))
bvr_points <- st_as_sf(
  bvr_points,
  coords = c("x", "y"),   # change to your column names
  crs = st_crs(thermal_masked_vec)    # match raster CRS!
)
# ccr
ccr_points <- data.frame(x = c(592750, 593050), y = c(4137200, 4139000))
ccr_points <- st_as_sf(
  ccr_points,
  coords = c("x", "y"),   # change to your column names
  crs = st_crs(thermal_masked_vec)    # match raster CRS!
)

################################################################################
# function to extract values and write out csv
################################################################################
get_vals <- function(points, thermal_data){
  vals <- st_extract(thermal_data["thermal_C"], points)
  vals_df <- data.frame(vals)
  if(dim(points)[1] > 1){
    vals_df <- vals_df %>%
      group_by(time) %>%
      summarize(mean_thermal_C = mean(thermal_C, na.rm = T))
    return(vals_df)
  } else {
    return(vals_df)
  }
}

################################################################################
# call functions and export csv
################################################################################
# set date of interest
#start_date <- "2013-04-19T00:00:00Z" # start of LS8
start_date <- "2025-01-19T00:00:00Z"
end_date <- "2025-11-05T00:00:00Z"
# call functions
data_ccr <- get_lst(ccr_box, ccr_box_utm, start_date, end_date)
thermal_masked_ccr <- water_mask(data_ccr)

# plot
ggplot() +
  geom_stars(data = thermal_masked_ccr["thermal_C"]) +
  facet_wrap(~time) +
  annotate("point", x = 592750, y = 4137200, color = 'red') +
  annotate("point", x = 593050, y = 4139000, color = 'red') +
  theme_classic() +
  scale_fill_viridis()
ccr_vals <- get_vals(ccr_points, thermal_masked_ccr)
write_csv(ccr_vals, "thermal_LS_CCR.csv")

