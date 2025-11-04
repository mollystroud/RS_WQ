################################################################################
# Code started by Molly Stroud on 10/20/25
################################################################################
# load in packages
require(pacman)
p_load('ggplot2', 'tidyverse', 'dplyr', 'patchwork')
################################################################################
# the below code is designed to explore correlations between filtered chla data,
# EXO data, and FLORA data
################################################################################

################################################################################
# EXO chla
################################################################################
# CCR
exo_ccr <- read_csv("ccre-waterquality_2021_2024.csv")
exo_ccr <- exo_ccr[!is.na(exo_ccr$EXOChla_ugL_1),] %>%
  select(DateTime, EXOChla_ugL_1) 
# average over each date
exo_ccr$DateTime <- as.Date(exo_ccr$DateTime)
exo_ccr <- exo_ccr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1)) %>%
  mutate(Latitude = 37.3697, Longitude =-79.958)
exo_ccr <- exo_ccr[!is.na(exo_ccr$mean_chla),]
exo_ccr$Reservoir <- "CCR"
exo_ccr_dates <- unique(as.Date(exo_ccr$DateTime))
# FCR
exo_fcr <- read_csv("fcre-waterquality_2018_2024.csv")
exo_fcr <- exo_fcr[!is.na(exo_fcr$EXOChla_ugL_1),] %>%
  select(DateTime, EXOChla_ugL_1) 
# average over each date
exo_fcr$DateTime <- as.Date(exo_fcr$DateTime)
exo_fcr <- exo_fcr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1)) %>%
  mutate(Latitude = 37.30325, Longitude =-79.8373)
exo_fcr <- exo_fcr[!is.na(exo_fcr$mean_chla),]
exo_fcr$Reservoir <- "FCR"
exo_fcr_dates <- unique(as.Date(exo_fcr$DateTime))
# BVR
exo_bvr <- read_csv("bvre-waterquality_2020_2024.csv")
exo_bvr <- exo_bvr[!is.na(exo_bvr$EXOChla_ugL_1.5),] %>%
  select(DateTime, EXOChla_ugL_1.5) 
# average over each date
exo_bvr$DateTime <- as.Date(exo_bvr$DateTime)
exo_bvr <- exo_bvr %>%
  group_by(DateTime) %>%
  summarize(mean_chla = mean(EXOChla_ugL_1.5)) %>%
  mutate(Latitude = 37.31288, Longitude =-79.8159)
exo_bvr <- exo_bvr[!is.na(exo_bvr$mean_chla),]
exo_bvr$Reservoir <- "BVR"
exo_bvr_dates <- unique(as.Date(exo_bvr$DateTime))
################################################################################

################################################################################
# filtered chla
################################################################################
chla <- read_csv("filt-chla_2014_2024.csv")
# new filtered chla 2025
chla_2025 <- read_csv("https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataNotYetUploadedToEDI/Raw_chla/Filt_chla_L1.csv")
chla <- data.frame(rbind(chla, chla_2025))
# only surface measurements
#chla <- chla[chla$Depth_m <= 0.1,]
# match with lat/long
sitelocs <- read_csv("site_descriptions.csv")
chla <- left_join(chla, sitelocs, by = c("Reservoir" = "Reservoir", "Site" = "Site"))

# only ccr, bvr, fcr
filtered_chla_ccr <- chla[chla$Reservoir == "CCR",] %>%
  filter(Site == 50)
filtered_chla_ccr$DateTime <- as.Date(filtered_chla_ccr$DateTime)
filtered_chla_fcr <- chla[chla$Reservoir == "FCR",] %>%
  filter(Site == 50)
filtered_chla_fcr$DateTime <- as.Date(filtered_chla_fcr$DateTime)
filtered_chla_bvr <- chla[chla$Reservoir == "BVR",] %>%
  filter(Site == 50)
filtered_chla_bvr$DateTime <- as.Date(filtered_chla_bvr$DateTime)

# within filtered, see how 0.1 and 1.6 match up
#filtered_chla_fcr <- filtered_chla_fcr %>%
  #group_by(DateTime, Depth_m) %>%
  #summarize(mean_chla = mean(Chla_ugL))

test <- filtered_chla_fcr[filtered_chla_fcr$Depth_m < 2,] %>%
  select(DateTime, Depth_m, mean_chla) %>%
  pivot_wider(names_from = Depth_m, values_from = mean_chla)

bvr_comp <- ggplot() +
  geom_point(data = test, 
             aes(x = `0.1`, y = `1.6`, color = DateTime)) +
  theme_classic() +
  labs(x = "0.1 m", y = '1.6 m', title = 'FCR') +
  geom_abline(slope = 1, intercept = 0)
bvr_comp
################################################################################

################################################################################
# plot EXO vs. filtered chla at 1.6 m
matched_ccr <- inner_join(filtered_chla_ccr, exo_ccr, by = "DateTime")
matched_fcr <- inner_join(filtered_chla_fcr, exo_fcr, by = "DateTime")
matched_bvr <- inner_join(filtered_chla_bvr, exo_bvr, by = "DateTime") #none

bvr_matchup <- ggplot(matched_bvr[matched_bvr$Depth_m == 1.6,], 
                      aes(x = Chla_ugL, y = mean_chla, color = DateTime)) +
  geom_point() + 
  theme_classic() +
  labs(x = 'Filtered Chla', y = 'EXO Chla', title = 'BVR, 1.6 m') +
  ylim(0, 90)
bvr_matchup

fcr_matchup <- ggplot(matched_fcr[matched_fcr$Depth_m == 1.6,], 
                      aes(x = Chla_ugL, y = mean_chla, color = DateTime)) +
  geom_point() + 
  theme_classic() +
  labs(x = 'Filtered Chla', y = 'EXO Chla', title = 'FCR, 1.6 m') +
  ylim(0, 100)
fcr_matchup
bvr_matchup + fcr_matchup

ggplot() +
  geom_point(data = matched_fcr[matched_fcr$Depth_m == 1.6,], 
             aes(x = Chla_ugL, y = mean_chla, color = 'darkblue'), alpha = 0.6) +
  geom_point(data = matched_bvr[matched_bvr$Depth_m == 1.6,],
             aes(x = Chla_ugL, y = mean_chla, color = 'red'), alpha = 0.6) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = 'Filtered Chla', y = 'EXO', color = element_blank()) +
  scale_color_manual(values =c('darkblue'='darkblue','red'='red'), labels = c('FCR','BVR')) +
  xlim(0, 100) + ylim(0, 100)

################################################################################
# subtract EXO from filtered
ggplot() +
  geom_point(data = matched_fcr[matched_fcr$Depth_m == 1.6,], 
             aes(x = DateTime, y = (Chla_ugL - mean_chla), color = 'darkblue'), alpha = 0.6) +
  geom_point(data = matched_bvr[matched_bvr$Depth_m == 1.6,],
             aes(x = DateTime, y = (Chla_ugL - mean_chla), color = 'red'), alpha = 0.6) +
  theme_classic() +
  geom_abline(slope = 0, intercept = 0) +
  labs(x = element_blank(), y = 'Filtered Chla - EXO', color = element_blank()) +
  scale_color_manual(values =c('darkblue'='darkblue','red'='red'), labels = c('FCR','BVR'))



################################################################################
# now add in flora
################################################################################
flora <- read_csv("fluoroprobe_2014_2024.csv")
flora$DateTime <- as.Date(flora$DateTime)
exo <- data.frame(rbind(exo_bvr, exo_ccr, exo_fcr))

# join FLORA and EXO
flora_exo <- inner_join(flora, exo, by = c("DateTime", "Reservoir"))
flora_exo <- flora_exo[flora_exo$Depth_m < 1.6 & flora_exo$Depth_m > 1.4,]
# plot
ggplot(flora_exo, aes(x = (TotalConc_ugL - Bluegreens_ugL), y = mean_chla, color = DateTime)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(y = "EXO Chla") +
  xlim(0, 65) + ylim(0, 65)

# looked at filtered chla
filtered_chla <- chla
filtered_chla$DateTime <- as.Date(filtered_chla$DateTime)
# 'surface'
flora_surface <- flora[flora$Depth_m > 0.48 & flora$Depth_m < 0.52,]
# match by date and reservoir
flora_filtered <- inner_join(flora_surface, filtered_chla, by = c("DateTime", "Reservoir"))

# only 'good' data
flora_filtered <- flora_filtered[flora_filtered$Flag_Chla_ugL == 0,]

# plot
ggplot(flora_filtered, aes(x = TotalConc_ugL, y = Chla_ugL, color = DateTime)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(y = "Filtered Chla")





