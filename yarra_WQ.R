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
#yarra_chla <- yarra_wq[yarra_wq$variable == "turbidity",]

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
yarra_box_small <- c(
  xmin = 145.895, 
  ymin = -37.680, 
  xmax = 145.905, 
  ymax = -37.670
)
# convert the bounding boxes to the correct UTM projection
yarra_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = 145.886, 
                                               ymin = -37.723, 
                                               xmax = 145.983, 
                                               ymax = -37.642), 
                                             crs = "EPSG:4326")), "EPSG:32655"))
yarra_box_utm_small <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(  xmin = 145.895, 
                                                 ymin = -37.680, 
                                                 xmax = 145.905, 
                                                 ymax = -37.670), 
                                             crs = "EPSG:4326")), "EPSG:32655"))

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
    ext_query("eo:cloud_cover" < 20) %>% # filter for cloud cover
    post_request() %>%
    items_fetch()
  dates <- sapply(items$features, function(x) x$properties$datetime)
}

# call function to get dates
yarra_hls_datetime <- get_dates(yarra_box_small, rev(yarra_dates)) # need to order in ascending chrono, rev if needed
yarra_hls_datetime

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
    max_dist = 1, # day difference
    distance_col = "day_diff"
  )
}

yarra_joined <- fuzzy_join_dates(yarra_hls_datetime, (yarra_dates))

# now join back with full chla data
yarra_chla$datetime <- as.Date(yarra_chla$datetime)
yarra_chla_joined <- yarra_chla[yarra_chla$datetime %in% yarra_joined$chla_dates,]
yarra_chla_joined <- left_join(yarra_chla_joined, yarra_joined, 
                             by = c("datetime" = "chla_dates"))
yarra_chla_joined <- unique(yarra_chla_joined)
yarra_chla_joined
################################################################################
# function to perform stac_search on specific dates
################################################################################
get_items <- function(joined_df, bbox){
  # grab items ON dates of interest
  dates <- unique(joined_df$HLS_dates)
  # get list of imagery
  items_list <- lapply(dates, function(d) {
    stac_search(q = s,
                collections = HLS_col,
                bbox = bbox,
                datetime = paste0(d,"T00:00:00Z", "/", d, "T23:59:59Z"),
                limit = 500) %>% 
      ext_query("eo:cloud_cover" < 20) %>%
      post_request()
  })
  
  # merge list
  items <- items_list[[1]]
  for(i in 2:length(items_list)){
    items$features <- c(items$features, items_list[[i]]$features)
  }
  
  # update metadata
  items$numberMatched <- length(items$features)
  items$numberReturned <- length(items$features)
  class(items) <- c("stac_item_collection", "list")
  # manually add projection info into each STAC item before passing to gdalcubes
  # for some reason HLS data doesn't have this automatically included?
  for (i in seq_along(items$features)) {
    items$features[[i]]$properties$`proj:epsg` <- 32655  # force correct projection
  }
  
  return(items)
}

items_yarra <- get_items(yarra_chla_joined, yarra_box_small)
items_yarra
################################################################################
# function to extract values HLSS or HLSL data based on lat/longs
################################################################################
get_vals <- function(HLStype, itemlist, joined_df, bbox_utm){
  # first define cube space
  # start/end based on dataframe
  start_date <- paste0(joined_df$HLS_dates[1], "T00:00:00Z")
  end_date <- paste0(joined_df[length(joined_df$HLS_dates),]$HLS_dates, "T00:00:00Z")
  # define the cube space
  cube <- cube_view(srs ="EPSG:32655",
                    extent = list(t0 = start_date, 
                                  t1 = end_date,
                                  left = bbox_utm[1], 
                                  right = bbox_utm[3],
                                  top = bbox_utm[4], 
                                  bottom = bbox_utm[2]),
                    dx = 30, # 30 m resolution
                    dy = 30, 
                    dt = 'P1D',
                    aggregation = "median", 
                    resampling = "average")
  # get only HLSL or HLSS stac image collection
  if(HLStype == "HLSS"){
    hls_items <- rev(itemlist$features)[sapply(rev(itemlist$features), function(f) f$collection == "HLSS30_2.0")]
    # sentinel collection
    col_S30 <- stac_image_collection(
      hls_items,
      asset_names = c("B02", "B03", "B04", "B8A"),
      url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
    )
    data <- raster_cube(image_collection = col_S30, 
                        view = cube)
    data <- rename_bands(data, B02 = "blue", B03 = "green", B04 = "red",
                         B8A = "NIR")
  } else {
    hls_items <- rev(itemlist$features)[sapply(rev(itemlist$features), function(f) f$collection == "HLSL30_2.0")]
    # landsat collection
    col_L30 <- stac_image_collection(
      hls_items,
      asset_names = c("B02", "B03", "B04", "B05"),
      url_fun = function(url) paste0("/vsicurl/", url)  # helps GDAL access
    )
    data <- raster_cube(image_collection = col_L30, 
                        view = cube)
    data <- rename_bands(data, B02 = "blue", B03 = "green", B04 = "red",
                         B05 = "NIR")
  }
  # get unique lat/longs
  points <- data.frame(unique(cbind(joined_df$Longitude, joined_df$Latitude)))
  points <- mutate(points, FID = row_number()) # add ID for merging later
  projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # project
  # make sf for extract_geom
  points_crs <- st_as_sf(x = points,                         
                         coords = c("X1", "X2"),
                         crs = projcrs)
  # extract band values at geometry points
  band_vals <- extract_geom(data, points_crs)
  band_vals <- left_join(band_vals, points, by = "FID") # add back lat/long
  band_vals$time <- as.Date(band_vals$time) # make Date
  
  return(band_vals)
}


# add lat/long
yarra_chla_joined$Latitude <- -37.67628
yarra_chla_joined$Longitude <- 145.90053
#yarra_chla_joined$Latitude[1] <- -37.67631
#yarra_chla_joined$Longitude[1] <- 145.90042
# dates need to be ordered correctly
yarra_chla_joined <- yarra_chla_joined[nrow(yarra_chla_joined):1, ]


# call function
yarra_vals_HLSS <- get_vals("HLSS", items_yarra, yarra_chla_joined, yarra_box_utm_small)
yarra_vals_HLSL <- get_vals("HLSL", items_yarra, yarra_chla_joined, yarra_box_utm_small)
yarra_vals <- data.frame(rbind(yarra_vals_HLSS, yarra_vals_HLSL)) # join
yarra_vals$time <- as.Date(yarra_vals$time)
yarra_alldata <- left_join(yarra_chla_joined, yarra_vals, # join with in situ
                         by = c("HLS_dates" = "time", "Latitude" = "X2", "Longitude" = "X1"))
yarra_alldata <- yarra_alldata[!is.na(yarra_alldata$FID),] # remove nas
yarra_alldata <- yarra_alldata[!duplicated(yarra_alldata), ] # remove duplicates
write_csv(yarra_alldata, "yarra_chla_matchups_2day.csv")


# CCR
yarra_alldata <- read_csv("yarra_chla_matchups_2day.csv") # read in data

yarra_alldata <- yarra_alldata[yarra_alldata$NIR > 0 & yarra_alldata$red > 0 &
                                 yarra_alldata$green > 0 & yarra_alldata$blue > 0,]
# create model 
model_yarra <- lm(log10(observation) ~ blue + green + red + NIR, data = yarra_alldata)
model_yarra_preds <- data.frame(cbind((predict(model_yarra)), yarra_alldata$observation))
summary(model_yarra) # summary stats
sqrt(mean(model_yarra$residuals^2)) # rmse
# plot
yarra_plot <- ggplot(model_yarra_preds, aes(x = 10^X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Yarra Reservoir") 
yarra_plot




