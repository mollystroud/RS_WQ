################################################################################
# Code started by Molly Stroud on 9/29/25
################################################################################
# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf')
# likely unused: 'imager', 'lubridate','xts','dygraphs', 'leaflet', 'cubelyr'
################################################################################
## the below code is designed to pull HLS imagery over a specified area and 
# mask out any land
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

# set date of interest
start_date <- "2013-01-01T00:00:00Z"
end_date <- "2014-12-31T00:00:00Z"

# set GDAL config options for HTTP / auth:
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEJAR", "/tmp/cookies.txt")
gdalcubes_set_gdal_config("GDAL_HTTP_COOKIEFILE", "/tmp/cookies.txt")
# and possibly disable listing on open
gdalcubes_set_gdal_config("CPL_VSIL_CURL_USE_HEAD", "FALSE")
gdalcubes_set_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN", "YES")

# grab items within dates of interest
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = fcr_box,
              datetime = paste(start_date, end_date, sep="/"),
              limit = 1000) %>%
  ext_query("eo:cloud_cover" < 20) %>% #filter for cloud cover
  post_request() %>%
  items_fetch()
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
                                left = fcr_box_utm[1], 
                                right = fcr_box_utm[3],
                                top = fcr_box_utm[4], 
                                bottom = fcr_box_utm[2]),
                  dx = 30, # 30 m resolution
                  dy = 30, 
                  dt = "P1D",
                  aggregation = "median", 
                  resampling = "average")

# make raster cube
# HLSS
data_S30 <- raster_cube(image_collection = col_S30, 
                     view = cube)
data_S30 <- rename_bands(data_S30, B02 = "blue", B03 = "green", B04 = "red",
                         B8A = "NIR", B11 = "SWIR")
#HLSL
data_L30 <- raster_cube(image_collection = col_L30, 
                        view = cube)
data_L30 <- rename_bands(data_L30, B02 = "blue", B03 = "green", B04 = "red",
                         B05 = "NIR", B06 = "SWIR")
# vnir (for future, once mask is created)
# HLSS
vnir_S30 <- data_S30 %>%
  select_bands(bands = c("red", "green", "blue", "NIR"))
#HLSL
vnir_S30 <- data_L30 %>%
  select_bands(bands = c("red", "green", "blue", "NIR"))

# calculate mNDWI for water classification
# HLSS
mndwi_S30 <- data_S30 %>%
  select_bands(bands = c("green", "SWIR")) %>%
  apply_pixel(expr = "(green-SWIR)/(green+SWIR)", names = "mNDWI")
# HLSL
mndwi_L30 <- data_L30 %>%
  select_bands(bands = c("green", "SWIR")) %>%
  apply_pixel(expr = "(green-SWIR)/(green+SWIR)", names = "mNDWI")

# open raster cubes as stars objects and merge!
mndwi_stars_S30 <- st_as_stars(mndwi_S30)
mndwi_stars_L30 <- st_as_stars(mndwi_L30)
mndwi_stars <- c(mndwi_stars_S30, mndwi_stars_L30, along = 3)

# remove dates with no imagery
arr <- mndwi_stars[[1]] # extract raw array (x, y, time)
non_na_counts <- apply(arr, 3, function(slice) sum(!is.na(slice))) # count non-NA pixels for each time
valid_idx <- which(non_na_counts > 0) # indices of slices that have at least one real value
# build cleaned object by stacking only valid slices
slices <- lapply(valid_idx, function(i) mndwi_stars[,,, i, drop = FALSE])
clean_mndwi_stars <- do.call(c, c(slices, along = "time"))
# viz
tm_shape(shp = clean_mndwi_stars) + 
  tm_raster(col = "mNDWI", palette = 'viridis')

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
  slice <- clean_mndwi_stars[,,,i, drop = TRUE] # get date
  vals <- as.vector(slice[[1]]) # vectorize
  vals_scaled <- (vals - min(vals, na.rm = TRUE)) / (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
  threshold_i <- otsu(matrix(vals_scaled, ncol = 1)) # calculate Otsu thresh
  thresholds[i] <- threshold_i
  mask_array[,,i] <- slice[[1]] > threshold_i # append
  print(threshold_i)
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
# rename
names(binary_stars) <- "watermask"
# viz
tm_shape(shp = binary_stars) + 
  tm_raster(col = "watermask", palette = '-viridis')

# next up! grab + merge VNIR data, then apply mask
# vnir (for future, once mask is created)
# HLSS
vnir_S30 <- data_S30 %>%
  select_bands(bands = c("red", "green", "blue", "NIR"))
#HLSL
vnir_L30 <- data_L30 %>%
  select_bands(bands = c("red", "green", "blue", "NIR"))

# open raster cubes as stars objects and merge!
vnir_stars_S30 <- st_as_stars(vnir_S30)
vnir_stars_L30 <- st_as_stars(vnir_L30)
vnir_stars <- c(vnir_stars_S30, vnir_stars_L30, along = 3)
# again, remove empty dates
# buildx cleaned object by stacking only valid slices
slices <- lapply(valid_idx, function(i) vnir_stars[,,, i, drop = FALSE])
clean_vnir_stars <- do.call(c, c(slices, along = "time"))

# mask out land
masked <- clean_vnir_stars[binary_stars]
# scale for viz
masked_scaled <- masked %>%
  mutate(across(everything(), ~ . * 0.0001)) %>%
  mutate(across(everything(), ~ pmin(pmax(., 0), 1)))  # clip to be safe
write_mdim(masked_scaled, "reservoirs_vnir_data/fcr_2013_2014.nc")
#test <- read_stars("reservoirs_vnir_data/fcr_2013_20215.nc")
# viz
tm_shape(shp = masked_scaled) + 
  tm_rgb(r = "red", g = "green", b = "blue",
         col.scale = tm_scale_rgb(max_color_value = .08))


# apply chla algorithm from chla_algorithm.R
coords <- expand.grid(
  x = st_get_dimension_values(masked, "x"),
  y = st_get_dimension_values(masked, "y"),
  time = st_get_dimension_values(masked, "time")
)
df <- data.frame(
  coords,
  red   = as.vector(masked$red),
  green = as.vector(masked$green),
  blue  = as.vector(masked$blue),
  NIR   = as.vector(masked$NIR)
)
df <- df[df$red > 0 & df$green > 0 & df$blue > 0 & df$NIR > 0,]
df <- na.omit(df)


#df$chla_pred <- exp(predict(gam_fit, newdata=df[4:7])) - 0.01
df$chla_pred_lm <- 10^(predict(model_fcr, df[4:7]))
df <- df[df$time != "2015-02-20",]
chla_stars <- st_as_stars(df, dims = c("x", "y", "time"), values = "chla_pred")
fcr_point <- data.frame(x = 603033.1, y = 4129159)
fcr_point <- st_as_sf(
  fcr_point,
  coords = c("x", "y"),   # change to your column names
  crs = st_crs(chla_stars)    # match raster CRS!
)
##########
# get closest values to dock
# convert stars to long format
chla_long <- as.data.frame(chla_stars, xy = TRUE, long = TRUE)

# keep only chla
chla_long <- chla_long %>%
  dplyr::select(x, y, time, chla_pred_lm) %>%
  filter(!is.na(chla_pred_lm))

# convert back to sf
chla_sf <- st_as_sf(chla_long,
                    coords = c("x", "y"),
                    crs = st_crs(chla_stars))

# order time
times <- sort(unique(chla_sf$time))

# extract nearest non-NA per date
result <- do.call(rbind, lapply(times, function(tt) {
  slice <- chla_sf %>% filter(time == tt)
  
  # nearest pixel with a non-NA value
  idx <- st_nearest_feature(fcr_point, slice)
  
  chosen <- slice[idx, ]
  chosen$time <- as.Date(chosen$time)
  chosen
}))

# get nearest value per date
result_df <- result %>%
  st_drop_geometry() %>%
  dplyr::select(time, value = chla_pred_lm)
result_df
############


# viz
ggplot() +
  geom_stars(data = chla_stars, aes(fill = chla_pred_lm)) +
  facet_wrap(~time, ncol = 4) +
  scale_fill_viridis_c(
    option = "viridis",
    na.value = "white", #+ #,
    trans = 'log') + # make NA cells white instead of gray
    theme_void() #+
  #geom_point(data = fcr_point, aes(x, y), color = 'red', size = 3)


# get mean chl-a values
# conver to df
df <- as.data.frame(chla_stars, long = TRUE)
df <- df[df$chla_pred_lm < 300 | is.na(df$chla_pred_lm) == TRUE,]
# compute mean for each date
mean_by_date <- df %>%
  group_by(time) %>%
  summarise(mean_value = mean(chla_pred_lm, na.rm = TRUE))
#mean_by_date <- mean_by_date[1:24,]
# output
mean_by_date$dock_value <- result_df$value
mean_by_date

write_csv(mean_by_date, "reservoirs_vnir_data/fcr_2013_2014.csv")





# read all together
values_files <- list.files('reservoirs_vnir_data/', pattern = '.csv', full.names = TRUE)
alldata <- data.frame()
for(file in values_files){
  data <- read_csv(file)
  alldata <- rbind(alldata, data)
}
alldata
alldata <- alldata[alldata$mean_value < 40 & alldata$dock_value < 40,]
write_csv(alldata, 'reservoirs_vnir_data/all_chla_est.csv')

ggplot() +
  geom_point(data = alldata, aes(x = time, y = mean_value), color = 'darkblue') +
  geom_line(data = alldata, aes(x = time, y = mean_value), color = 'darkblue') +
  geom_point(data = alldata, aes(x = time, y = dock_value), color = 'red') +
  geom_line(data = alldata, aes(x = time, y = dock_value), color = 'red') +
  theme_classic()



