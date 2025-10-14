# Observing lakes with SWOT

################################################################################
# libraries
################################################################################
library(tidyverse)
library(sf)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(viridis)
################################################################################
# plot LakeSP product over roanoke reservoirs
################################################################################

################################################################################
# read in all SWOT data and filter
################################################################################
lakesp_files <- list.files("/Users/mollystroud/Desktop/postdoc/SWOT", pattern = '.shp$', 
                           full.names = T)
lakesp_files
# read in each file and append to one big dataframe
lakedata <- data.frame()
for(file in lakesp_files[4]){
  print(file)
  data <- read_sf(file)
  sf <- data
  #data <- data[data$wse > 1,]
  #data <- select(data, -one_of("qual_f_b"))
  #lakedata <- rbind(lakedata, data)
}

#plot(sf)


# ccr
ccr <- sf[sf$lake_name == 'CARVINS COVE RESERVOIR',]
# get outline
ccr_outline <- read_sf("/Users/mollystroud/Desktop/postdoc/SWOT/Reservoir_(2024_Final_WQA_IR_Assessment)/Reservoir_(2024_Final_WQA_IR_Assessment).shp")

ccr_map <- ggplot() + 
  geom_sf(data = ccr, aes(fill = wse), color = NA) + 
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Carvin Cove Reservoir",], 
          color = 'black', fill = NA, size = 0.5) +
  coord_sf() + 
  scale_fill_viridis() +
  theme_void() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = 'WSE (m)') +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")
ccr_map


st_write(ccr, "/Users/mollystroud/Desktop/postdoc/SWOT/ccr_map.shp")
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/ccr_map.pdf", plot = ccr_map, width = 5, height = 5, units = c("in"))

# fcr
fcr <- sf[sf$lake_name == 'FALLING CREEK RESERVOIR',]
fcr_map <- ggplot() + 
  geom_sf(data = fcr, aes(fill = wse), color = NA) + 
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Falling Creek Reservoir" & 
                               ccr_outline$AU_Categor == "3A",],
          color = 'black', fill = NA, size = 0.5) + 
  coord_sf() +
  scale_fill_viridis() +
  theme_void() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = 'WSE (m)') +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")
fcr_map
st_write(fcr, "/Users/mollystroud/Desktop/postdoc/SWOT/fcr_map.shp")
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/fcr_map.pdf", plot = fcr_map, width = 5, height = 5, units = c("in"))
# bdr
bdr <- sf[sf$lake_id == '7320364932',]
bdr_map <- ggplot() + 
  geom_sf(data = bdr, aes(fill = wse), color = NA) + 
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Beaverdam Reservoir (XKD)",], 
          color = 'black', fill = NA, size = 0.5) +
  coord_sf() +
  scale_fill_viridis() +
  theme_void() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = 'WSE (m)') +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")
bdr_map

st_write(bdr, "/Users/mollystroud/Desktop/postdoc/SWOT/bdr_map.shp")
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/bdr_map.pdf", plot = bdr_map, width = 5, height = 5, units = c("in"))

ccr_map + fcr_map + bdr_map












library(ncdf4)
library(raster)
library(sf)
library(terra)
library(ggplot2)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(elevatr)
library(terra)
library(sf)
library(ggnewscale)
library(stars)
library(usmap)


################################################################################ 
# Open PIXC
################################################################################
test <- nc_open('/Users/mollystroud/Desktop/postdoc/SWOT/SWOT_L2_HR_PIXC_028_397_220L_20250217T180033_20250217T180044_PIC2_01.nc')

#names(test$var)
# get vars
lon <- ncvar_get(test, "pixel_cloud/longitude")
lat <- ncvar_get(test, "pixel_cloud/latitude")
wse <- ncvar_get(test, "pixel_cloud/height")
class <- ncvar_get(test, "pixel_cloud/classification")
qual <- ncvar_get(test, "pixel_cloud/classification_qual")
# close
nc_close(test)
# merge
df <- data.frame(cbind(lon, lat, wse, class, qual))
# clean
df <- df[df$class == 4,]
#df <- df[df$class > 2,] # for all water types
df <- df[df$qual <= 2,]
#df <- df[df$wse < 800,]
#df <- df[df$wse > 400,]

################################################################################ 
# Plot PIXC 
################################################################################ 
# plot CCR (SWOT_L2_HR_PIXC_028_397_220L_20250217T180033_20250217T180044_PIC2_01.nc)
ccr_df <- df[df$lon > -79.981728 & df$lon < -79.942552 &
               df$lat > 37.367522 & df$lat < 37.407255,]

swot_ccr <- ggplot() +
  geom_point(data = ccr_df, aes(x = lon, y = lat, color = wse), size = 0.1) +
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Carvin Cove Reservoir",], 
          color = 'black', fill = NA, size = 0.5) +
  scale_color_viridis() +
  labs(x = element_blank(), y = element_blank(), 
       color = 'WSE (m)') +
  coord_sf(crs = 4326) +
  theme_void() +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")

swot_ccr
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_ccr_map.pdf", plot = swot_ccr, width = 5, height = 5, units = c("in"))
# make ccr raster
# set up a raster geometry, here deriving an extent from your data
e <- ext(apply(ccr_df[,1:2], 2, range))
# set up the raster, for example
r <- rast(e, ncol=130, nrow=130)
# you need to provide a function 'fun' for when there are multiple points per cell
ccr_rast <- rasterize(ccr_df[, 1:2], r, ccr_df[,3], fun=mean)
plot(x)
writeRaster(ccr_rast, "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_ccr_map.tif", overwrite = T)

test <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_ccr_map.tif")
plot(test)


# plot FCR (SWOT_L2_HR_PIXC_028_397_220L_20250217T180033_20250217T180044_PIC2_01.nc)
fcr_df <- df[df$lon > -79.840037 & df$lon < -79.833651 &
               df$lat > 37.301435 & df$lat < 37.311487,]
swot_fcr <- ggplot() +
  geom_point(data = fcr_df, aes(x = lon, y = lat, color = wse), size = 1) +
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Falling Creek Reservoir" & 
                               ccr_outline$AU_Categor == "3A",],
          color = 'black', fill = NA, size = 0.5) + 
  scale_color_viridis() +
  labs(x = element_blank(), y = element_blank(), 
       color = 'WSE (m)') +
  coord_sf(crs = 4326) +
  theme_void() +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")
swot_fcr
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_fcr_map.pdf", plot = swot_fcr, width = 5, height = 5, units = c("in"))
# make fcr raster
# set up a raster geometry, here deriving an extent from your data
e <- ext(apply(fcr_df[,1:2], 2, range))
# set up the raster, for example
r <- rast(e, ncol=20, nrow=20)
# you need to provide a function 'fun' for when there are multiple points per cell
fcr_rast <- rasterize(fcr_df[, 1:2], r, fcr_df[,3], fun=mean)
plot(fcr_rast)
writeRaster(fcr_rast, "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_fcr_map.tif", overwrite = T)

test <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_fcr_map.tif")
plot(test)


# BDR 
bdr_df <- df[df$lon > -79.827088 & df$lon < -79.811865 &
               df$lat > 37.311798 & df$lat < 37.321694,]
swot_bdr <- ggplot() +
  geom_point(data = bdr_df, aes(x = lon, y = lat, color = wse), size = 1) +
  geom_sf(data = ccr_outline[ccr_outline$WaterName == "Beaverdam Reservoir (XKD)",], 
          color = 'black', fill = NA, size = 0.5) +
  scale_color_viridis() +
  labs(x = element_blank(), y = element_blank(), 
       color = 'WSE (m)') +
  coord_sf(crs = 4326) +
  theme_void() +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = "black")
swot_bdr
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_bdr_map.pdf", plot = swot_bdr, width = 5, height = 5, units = c("in"))

# make bdr raster
# set up a raster geometry, here deriving an extent from your data
e <- ext(apply(bdr_df[,1:2], 2, range))
# set up the raster, for example
r <- rast(e, ncol=30, nrow=30)
# you need to provide a function 'fun' for when there are multiple points per cell
bdr_rast <- rasterize(bdr_df[, 1:2], r, bdr_df[,3], fun=mean)
plot(bdr_rast)
writeRaster(bdr_rast, "/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_bdr_map.tif", overwrite = T)

test <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/swot_pixc_bdr_map.tif")
plot(test)

swot_ccr + swot_fcr + swot_bdr


# thermal data
# ccr
ccr_thermal <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/ccr_thermal_02_2025.tif")
ccr_outline_crs <- st_transform(ccr_outline, crs(ccr_thermal))
ccr_thermal_df <- as.data.frame(ccr_thermal, xy = TRUE)
ccr_thermal_plot <- ggplot() +
  geom_raster(data = ccr_thermal_df, aes(x = x, y = y, fill = (ST_B10 * 0.00341802 - 124.15))) +
  scale_fill_viridis() +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Carvin Cove Reservoir",], 
          color = 'black', fill = NA, size = 0.5) +
  theme_void() +
  labs(fill = "Temp (Celcius)")
ccr_thermal_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/thermal_ccr.pdf",
       plot = ccr_thermal_plot, width = 5, height = 5, units = c("in"))

# fcr + bvr
fcr_thermal <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/fcr_bdr_thermal_04_2025.tif")
#ccr_outline_crs <- st_transform(ccr_outline, crs(ccr_thermal))
fcr_thermal_df <- as.data.frame(fcr_thermal, xy = TRUE)
fcr_thermal_plot <- ggplot() +
  geom_raster(data = fcr_thermal_df, aes(x = x, y = y, fill = (ST_B10 * 0.00341802 - 124.15))) +
  scale_fill_viridis() +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Beaverdam Reservoir (XKD)",], 
          color = 'black', fill = NA, size = 0.5) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Falling Creek Reservoir" & 
                                   ccr_outline_crs$AU_Categor == "3A",],
          color = 'black', fill = NA, size = 0.5) + 
  
  theme_void() +
  labs(fill = "Temp (Celcius)")
fcr_thermal_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/thermal_fcr_bvr.pdf",
       plot = fcr_thermal_plot, width = 5, height = 5, units = c("in"))



# hls data
# ccr
ccr_ndwi <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/ccr_ndwi.tif")
ccr_outline_crs <- st_transform(ccr_outline, crs(ccr_ndwi))
ccr_ndwi_df <- as.data.frame(ccr_ndwi, xy = TRUE)
ccr_ndwi_plot <- ggplot() +
  geom_raster(data = ccr_ndwi_df, aes(x = x, y = y, fill = NDWI)) +
  scale_fill_gradient2(low = 'white', mid = '#0096C7', high = '#023E8A', limits = c(-1, 1)) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Carvin Cove Reservoir",], 
          color = 'black', fill = NA, size = 0.5) +
  theme_void()
ccr_ndwi_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/ndwi_ccr.pdf",
       plot = ccr_ndwi_plot, width = 5, height = 5, units = c("in"))

# fcr + bvr
fcr_bdr_ndwi <- raster("/Users/mollystroud/Desktop/postdoc/SWOT/fcr_bdr_ndwi.tif")
fcr_outline_crs <- st_transform(ccr_outline, crs(fcr_bdr_ndwi))
fcr_bdr_ndwi_df <- as.data.frame(fcr_bdr_ndwi, xy = TRUE)
fcr_ndwi_plot <- ggplot() +
  geom_raster(data = fcr_bdr_ndwi_df, aes(x = x, y = y, fill = NDWI)) +
  scale_fill_gradient2(low = 'white', mid = '#0096C7', high = '#023E8A', limits = c(-1, 1)) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Beaverdam Reservoir (XKD)",], 
          color = 'black', fill = NA, size = 0.5) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Falling Creek Reservoir" & 
                                   ccr_outline_crs$AU_Categor == "3A",],
          color = 'black', fill = NA, size = 0.5) + 
  theme_void()
fcr_ndwi_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/ndwi_fcr_bvr.pdf",
       plot = fcr_ndwi_plot, width = 5, height = 5, units = c("in"))



# planet data
library(terrainr)
# ccr
ccr_planet <- stack("/Users/mollystroud/Desktop/postdoc/SWOT/planet_ccr_2025_10_05.tif")
ccr_outline_crs <- st_transform(ccr_outline, crs(ccr_planet))
ccr_planet_df <- as.data.frame(ccr_planet, xy = TRUE)
colnames(ccr_planet_df) <- c("x", "y", "red", "green", "blue", "nir")
ccr_planet_plot <- ggplot() +
  geom_spatial_rgb(
    data = ccr_planet_df,
    mapping = aes(
      x = x,
      y = y,
      r = red,
      g = green,
      b = blue)) +
  #scale_fill_gradient2(low = 'white', mid = '#0096C7', high = '#023E8A', limits = c(-1, 1)) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Carvin Cove Reservoir",], 
          color = 'white', fill = NA, size = 0.5) +
  theme_void()
ccr_planet_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/planet_ccr.pdf",
       plot = ccr_planet_plot, width = 5, height = 5, units = c("in"))

# fcr + bvr
fcr_bdr_planet <- stack("/Users/mollystroud/Desktop/postdoc/SWOT/planet_fcr_bdr_2025_10_05.tif")
#fcr_outline_crs <- st_transform(ccr_outline, crs(fcr_bdr_ndwi))
fcr_bdr_planet_df <- as.data.frame(fcr_bdr_planet, xy = TRUE)
colnames(fcr_bdr_planet_df) <- c("x", "y", "red", "green", "blue", "nir")

fcr_planet_plot <- ggplot() +
  geom_spatial_rgb(
    data = fcr_bdr_planet_df,
    mapping = aes(
      x = x,
      y = y,
      r = red,
      g = green,
      b = blue)) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Beaverdam Reservoir (XKD)",], 
          color = 'white', fill = NA, size = 0.5) +
  geom_sf(data = ccr_outline_crs[ccr_outline_crs$WaterName == "Falling Creek Reservoir" & 
                                   ccr_outline_crs$AU_Categor == "3A",],
          color = 'white', fill = NA, size = 0.5) + 
  theme_void()
fcr_planet_plot
ggsave(filename = "/Users/mollystroud/Desktop/postdoc/SWOT/Figs for Illustrator/planet_fcr_bvr.pdf",
       plot = fcr_planet_plot, width = 5, height = 5, units = c("in"))



