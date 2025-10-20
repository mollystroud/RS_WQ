################################################################################
# Code started by Molly Stroud on 10/10/25
################################################################################

# load in packages
require(pacman)
p_load('earthdatalogin', 'rstac', 'terra', 'stars', 'ggplot2', 'tidyterra', 
       'viridis', 'EBImage', 'gdalcubes', 'tmap', 'dplyr', 'tidyverse', 'sf',
       'fuzzyjoin', 'patchwork')
################################################################################
# the below code is designed to match up in situ chl-a data over CCR, FCR, and
# BVR with HLS data
################################################################################
# read in filtered chl-a data from reservoirs
################################################################################
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
################################################################################
# or, EXO chla
################################################################################
# CCR
chla_ccr <- read_csv("ccre-waterquality_2021_2024.csv")
chla_ccr <- chla[!is.na(chla$EXOChla_ugL_1),] %>%
  select(DateTime, EXOChla_ugL_1) 
# average over each date
chla_ccr$DateTime <- as.Date(chla_ccr$DateTime)
chla_ccr <- chla_ccr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1)) %>%
  mutate(Latitude = 37.3697, Longitude =-79.958)
chla_ccr <- chla_ccr[!is.na(chla_ccr$mean_chla),]
ccr_dates <- unique(as.Date(chla_ccr$DateTime))
# FCR
chla_fcr <- read_csv("fcre-waterquality_2018_2024.csv")
chla_fcr <- chla_fcr[!is.na(chla_fcr$EXOChla_ugL_1),] %>%
  select(DateTime, EXOChla_ugL_1) 
# average over each date
chla_fcr$DateTime <- as.Date(chla_fcr$DateTime)
chla_fcr <- chla_fcr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1)) %>%
  mutate(Latitude = 37.30325, Longitude =-79.8373)
chla_fcr <- chla_fcr[!is.na(chla_fcr$mean_chla),]
fcr_dates <- unique(as.Date(chla_fcr$DateTime))
# BVR
chla_bvr <- read_csv("bvre-waterquality_2020_2024.csv")
chla_bvr <- chla_bvr[!is.na(chla_bvr$EXOChla_ugL_1.5),] %>%
  select(DateTime, EXOChla_ugL_1.5) 
# average over each date
chla_bvr$DateTime <- as.Date(chla_bvr$DateTime)
chla_bvr <- chla_bvr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1.5)) %>%
  mutate(Latitude = 37.31288, Longitude =-79.8159)
chla_bvr <- chla_bvr[!is.na(chla_bvr$mean_chla),]
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
ccr_hls_datetime <- get_dates(ccr_box, ccr_dates)
fcr_hls_datetime <- get_dates(fcr_box, fcr_dates)
bvr_hls_datetime <- get_dates(bvr_box, bvr_dates)

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

ccr_joined <- fuzzy_join_dates(ccr_hls_datetime, ccr_dates)
fcr_joined <- fuzzy_join_dates(fcr_hls_datetime, fcr_dates)
bvr_joined <- fuzzy_join_dates(bvr_hls_datetime, bvr_dates)
# for EXO data
ccr_joined <- ccr_joined[ccr_joined$chla_dates == ccr_joined$HLS_dates,]
fcr_joined <- fcr_joined[fcr_joined$chla_dates == fcr_joined$HLS_dates,]
bvr_joined <- bvr_joined[bvr_joined$chla_dates == bvr_joined$HLS_dates,]



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
    items$features[[i]]$properties$`proj:epsg` <- 32617  # force correct projection
  }
  
  return(items)
}

items_ccr <- get_items(ccr_chla_joined, ccr_box)
items_fcr <- get_items(fcr_chla_joined, fcr_box)
items_bvr <- get_items(bvr_chla_joined, bvr_box)

################################################################################
# function to extract values HLSS or HLSL data based on lat/longs
################################################################################
get_vals <- function(HLStype, itemlist, joined_df, bbox_utm){
  # first define cube space
  # start/end based on dataframe
  start_date <- paste0(joined_df$HLS_dates[1], "T00:00:00Z")
  end_date <- paste0(joined_df[length(joined_df$HLS_dates),]$HLS_dates, "T00:00:00Z")
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
                    dt = 'P1D',
                    aggregation = "median", 
                    resampling = "average")
  # get only HLSL or HLSS stac image collection
  if(HLStype == "HLSS"){
    hls_items <- itemlist$features[sapply(itemlist$features, function(f) f$collection == "HLSS30_2.0")]
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
    hls_items <- itemlist$features[sapply(itemlist$features, function(f) f$collection == "HLSL30_2.0")]
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

# ccr
ccr_vals_HLSS <- get_vals("HLSS", items_ccr, ccr_chla_joined, ccr_box_utm)
ccr_vals_HLSL <- get_vals("HLSL", items_ccr, ccr_chla_joined, ccr_box_utm)
ccr_vals <- data.frame(rbind(ccr_vals_HLSS, ccr_vals_HLSL)) # join
ccr_vals$time <- as.Date(ccr_vals$time)
ccr_alldata <- left_join(ccr_chla_joined, ccr_vals, # join with in situ
                         by = c("HLS_dates" = "time", "Latitude" = "X2", "Longitude" = "X1"))
ccr_alldata <- ccr_alldata[!is.na(ccr_alldata$FID),] # remove nas
ccr_alldata <- ccr_alldata[!duplicated(ccr_alldata), ] # remove duplicates
#write_csv(ccr_alldata, "EXO_chla_ccr_matchups.csv")

# fcr
fcr_vals_HLSS <- get_vals("HLSS", items_fcr, fcr_chla_joined, fcr_box_utm)
fcr_vals_HLSL <- get_vals("HLSL", items_fcr, fcr_chla_joined, fcr_box_utm)
fcr_vals <- rbind(fcr_vals_HLSS, fcr_vals_HLSL) # join
fcr_alldata <- left_join(fcr_chla_joined, fcr_vals, # join with in situ
                         by = c("HLS_dates" = "time", "Latitude" = "X2", "Longitude" = "X1"))
fcr_alldata <- fcr_alldata[!is.na(fcr_alldata$FID),]
fcr_alldata <- fcr_alldata[!duplicated(fcr_alldata), ]
#write_csv(fcr_alldata, "EXO_chla_fcr_matchups.csv")

# bvr
bvr_vals_HLSS <- get_vals("HLSS", items_bvr, bvr_chla_joined, bvr_box_utm)
bvr_vals_HLSL <- get_vals("HLSL", items_bvr, bvr_chla_joined, bvr_box_utm)
bvr_vals <- rbind(bvr_vals_HLSS, bvr_vals_HLSL) #join
bvr_alldata <- left_join(bvr_chla_joined, bvr_vals, # join with in situ
                         by = c("HLS_dates" = "time", "Latitude" = "X2", "Longitude" = "X1"))
bvr_alldata <- bvr_alldata[!is.na(bvr_alldata$FID),]
bvr_alldata <- bvr_alldata[!duplicated(bvr_alldata), ]
write_csv(bvr_alldata, "EXO_chla_bvr_matchups.csv")



# now model time
ccr_alldata <- read_csv("EXO_chla_ccr_matchups.csv")
# remove data with weird vals
ccr_alldata <- ccr_alldata[ccr_alldata$blue > 0 & ccr_alldata$green > 0 &
                             ccr_alldata$red > 0 & ccr_alldata$NIR > 0,]
# CCR

# LOOCV TEST
# Initialize variables to store RMSE values and actual vs predicted values
rmse <- numeric()
actual_vs_predicted <- data.frame(Actual = numeric(), Predicted = numeric())

# Perform Leave-One-Out Cross Validation
for (i in 1:nrow(bvr_alldata)) {
  # Exclude the ith row
  test_data <- bvr_alldata[i, ]
  train_data <- bvr_alldata[-i, ]
  
  # Train the model
  model <- lm(Chla_ugL ~ green + blue + red + NIR, data = train_data)
  
  # Make predictions on the test data
  predictions <- predict(model, newdata = test_data)
  
  # Calculate RMSE
  rmse[i] <- sqrt(mean((test_data$Chla_ugL - predictions)^2))
  
  # Store actual vs predicted values
  actual_vs_predicted <- rbind(actual_vs_predicted, data.frame(Actual = test_data$Chla_ugL, Predicted = predictions))
}

# Average RMSE across all folds
average_rmse <- round(mean(rmse), 2)
actual_vs_predicted %>%
  mutate(rmse = sqrt((Actual - Predicted)^2))


# GAM
library(mgcv)
ccr_alldata$log_chl <- log(ccr_alldata$mean_chla + 0.01)  # avoid zeros

# Helper function to fit GAM with specified k
fit_gam <- function(k_band) {
  gam(log_chl ~ 
        s(red,  k = k_band) +
        s(green,  k = k_band) +
        s(blue,  k = k_band) +
        s(NIR, k = k_band),
      data = ccr_alldata,
      method = "REML")
}


# Try a few complexity settings
gam1 <- fit_gam(k_band = 4)
gam2 <- fit_gam(k_band = 6)
gam3 <- fit_gam(k_band = 8)
gam4 <- fit_gam(k_band = 10)
gam5 <- fit_gam(k_band = 12)

# Compare using AIC
AIC(gam1, gam2, gam3, gam4, gam5)
# Check smoothness diagnostics
# gam.check(gam4)
# Summary of best model
summary(gam4)

ccr_alldata$pred_chl <- exp(predict(gam4, ccr_alldata)) - 0.01
# Performance metrics
rmse <- sqrt(mean((ccr_alldata$pred_chl - ccr_alldata$mean_chla)^2))
r2 <- cor(ccr_alldata$pred_chl, ccr_alldata$mean_chla)^2

cat("RMSE:", rmse, "\nR2:", r2, "\n")

par(mfrow = c(2,3))
plot(gam2, shade=TRUE)

# plot
ccr_plot <- ggplot(ccr_alldata, aes(x = pred_chl, y = mean_chla)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "EXO Chl-a estimates, CCR") #+
#annotate("text", x = 2, y = 8, label = "R2 = 0.76, RMSE = 1.14")
ccr_plot

# NOW TRAIN/TEST

ccr_alldata <- ccr_alldata[order(ccr_alldata$DateTime), ]  # sort by time

train_idx <- 1:floor(0.7 * nrow(ccr_alldata))
test_idx <- (floor(0.7 * nrow(ccr_alldata)) + 1):nrow(df)

train <- ccr_alldata[train_idx, ]
test  <- ccr_alldata[test_idx, ]

# Fit GAM
gam_fit <- gam(log_chl ~ s(red) + s(green) + s(blue) + s(NIR), data=train, method="REML")

# Predict on held-out set
test$pred <- exp(predict(gam_fit, newdata=test)) - 0.01

# Evaluate
rmse <- sqrt(mean((test$pred - test$mean_chla)^2))
r2   <- cor(test$pred, test$mean_chla^2)
cat("RMSE:", rmse, " R2:", r2, "\n")
ggplot(test, aes(x = pred, y = mean_chla)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0)


# FCR
model_fcr <- lm(Chla_ugL ~ blue + green + red + NIR, data = fcr_alldata)
model_fcr_preds <- data.frame(cbind(predict(model_fcr), fcr_alldata$Chla_ugL))
summary(model_fcr)
sqrt(mean(model_fcr$residuals^2))
# plot
fcr_plot <- ggplot(model_fcr_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Filtered Chl-a estimates, FCR") +
  annotate("text", x = 5, y = 25, label = "R2 = 0.77, RMSE = 5.15")
fcr_plot

# BVR
model_bvr <- lm(Chla_ugL ~ blue + green + red + NIR, data = bvr_alldata)
model_bvr_preds <- data.frame(cbind(predict(model_bvr), bvr_alldata$Chla_ugL))
summary(model_bvr)
sqrt(mean(model_bvr$residuals^2))
# plot
bvr_plot <- ggplot(model_bvr_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Filtered Chl-a estimates, BVR") +
  annotate("text", x = 2.5, y = 12.5, label = "R2 = 0.89, RMSE = 1.6")
bvr_plot

ccr_plot + fcr_plot + bvr_plot




#




