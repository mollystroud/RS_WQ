################################################################################
# Code started by Molly Stroud on 10/20/25
################################################################################
# load in packages
require(pacman)
p_load('mgcv', 'ggplot2', 'tidyverse', 'dplyr', 'patchwork')

################################################################################
# the below code is designed to use the match-ups from chla_data_match.R to
# create accurate water quality estimation algorithms
# and, to test how surface chla measurements compare to measurements at depth
################################################################################

# quick comparison of filtered and EXO to see relationship
chla <- read_csv("filt-chla_2014_2024.csv")
# new filtered chla 2025
chla_2025 <- read_csv("https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataNotYetUploadedToEDI/Raw_chla/Filt_chla_L1.csv")
chla <- data.frame(rbind(chla, chla_2025))
# only surface measurements
chla <- chla[chla$Depth_m <= 0.1,]
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

# within filtered
test <- filtered_chla_ccr %>%
  group_by(DateTime, Depth_m) %>%
  summarize(mean_chla = mean(Chla_ugL))
ccr_comp <- ggplot() +
  geom_line(data = test[#test$DateTime < "2019-01-01" & 
                          #test$DateTime > "2015-01-01" &
                          test$Depth_m <= 6,], 
             aes(x = DateTime, y = mean_chla, color = Depth_m, group = Depth_m),
            linewidth = 1) +
  theme_classic() +
  labs(x = element_blank(), y = 'Chl-a (ugL)', title = 'BVR')
ccr_comp
ccr_comp + fcr_comp + bvr_comp

# other dfs are from chla_data_match.R
matched_ccr <- inner_join(filtered_chla_ccr, chla_ccr, by = "DateTime")
matched_fcr <- inner_join(filtered_chla_fcr, chla_fcr, by = "DateTime")
matched_bvr <- inner_join(filtered_chla_bvr, chla_bvr, by = "DateTime") #none

ccr_matchup <- ggplot(matched_ccr, aes(x = Chla_ugL, y = mean_chla)) +
  geom_point() + 
  theme_classic() +
  labs(x = 'Filtered Chla', y = 'EXO Chla', titled = 'CCR')
fcr_matchup <- ggplot(matched_fcr, aes(x = Chla_ugL, y = mean_chla)) +
  geom_point() + 
  theme_classic() +
  labs(x = 'Filtered Chla', y = 'EXO Chla', title = 'FCR')
ccr_matchup + fcr_matchup

################################################################################
# create simple linear estimation algorithm with filtered chla
################################################################################
# CCR
ccr_alldata <- read_csv("filtered_chla_ccr_matchups_2day.csv") # read in data
ccr_alldata <- ccr_alldata[ccr_alldata$Flag_Chla_ugL < 2,]
# create model 
model_ccr <- lm((Chla_ugL) ~ blue + green + red + NIR, data = ccr_alldata)
model_ccr_preds <- data.frame(cbind((predict(model_ccr)), ccr_alldata$Chla_ugL))
summary(model_ccr) # summary stats
sqrt(mean(model_ccr$residuals^2)) # rmse
# plot
ccr_plot <- ggplot(model_ccr_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Filtered Chl-a estimates, CCR") 
ccr_plot

# FCR
fcr_alldata <- read_csv("filtered_chla_fcr_matchups.csv") # read in data
# remove really high #s
fcr_alldata <- fcr_alldata[fcr_alldata$NIR < 2000,]
# create model
model_fcr <- lm(log10(Chla_ugL) ~ blue + green + red + NIR, data = fcr_alldata)
model_fcr_preds <- data.frame(cbind(10^(predict(model_fcr)), fcr_alldata$Chla_ugL))
summary(model_fcr) # summary stats
sqrt(mean(model_fcr$residuals^2)) # rmse
# plot
fcr_plot <- ggplot(model_fcr_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Filtered Chl-a estimates, FCR") 
fcr_plot

# BVR
bvr_alldata <- read_csv("filtered_chla_bvr_matchups.csv") # read in data
# create model
model_bvr <- lm(log10(Chla_ugL) ~ blue + green + red + NIR, data = bvr_alldata)
model_bvr_preds <- data.frame(cbind(10^(predict(model_bvr)), (bvr_alldata$Chla_ugL)))
summary(model_bvr) # summary stats
sqrt(mean(model_bvr$residuals^2)) # rmse
# plot
bvr_plot <- ggplot(model_bvr_preds, aes(x = (X1), y = (X2))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "Filtered Chl-a estimates, BVR") #+
  annotate("text", x = 2.5, y = 12.5, label = "R2 = 0.89, RMSE = 1.6")
bvr_plot

ccr_plot + fcr_plot + bvr_plot


# LOOCV TEST if interested in seeing
# Initialize variables to store RMSE values and actual vs predicted values
rmse <- numeric()
actual_vs_predicted <- data.frame(Actual = numeric(), Predicted = numeric())

# Perform Leave-One-Out Cross Validation
for (i in 1:nrow(bvr_alldata)) {
  # Exclude the ith row
  test_data <- bvr_alldata[i, ]
  train_data <- bvr_alldata[-i, ]
  
  # Train the model
  model <- lm((Chla_ugL) ~ green + blue + red + NIR, data = train_data)
  
  # Make predictions on the test data
  predictions <- predict(model, newdata = test_data)
  
  # Calculate RMSE
  rmse[i] <- sqrt(mean((test_data$Chla_ugL - predictions)^2))
  
  # Store actual vs predicted values
  actual_vs_predicted <- rbind(actual_vs_predicted, 
                               data.frame(Actual = test_data$Chla_ugL, Predicted = predictions))
}

# Average RMSE across all folds
average_rmse <- round(mean(rmse), 2)
average_rmse
actual_vs_predicted %>%
  mutate(rmse = sqrt((Actual - Predicted)^2))

################################################################################
# test out using generalized additive model (GAM) for EXO data
################################################################################
# now model time
ccr_EXO <- read_csv("EXO_chla_ccr_matchups.csv")
fcr_EXO <- read_csv("EXO_chla_fcr_matchups.csv")
bvr_EXO <- read_csv("EXO_chla_bvr_matchups.csv")

# remove data with weird vals (negatives likely caused by clouds)
ccr_EXO <- ccr_EXO[ccr_EXO$blue > 0 & ccr_EXO$green > 0 &
                     ccr_EXO$red > 0 & ccr_EXO$NIR > 0,]

# GAM
ccr_EXO$log_chl <- log(ccr_EXO$mean_chla + 0.01)  # avoid zeros

# Helper function to fit GAM with specified k
fit_gam <- function(k_band) {
  gam(log_chl ~ 
        s(red,  k = k_band) +
        s(green,  k = k_band) +
        s(blue,  k = k_band) +
        s(NIR, k = k_band),
      data = ccr_EXO,
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

ccr_EXO$pred_chl <- exp(predict(gam4, ccr_EXO)) - 0.01
# Performance metrics
rmse <- sqrt(mean((ccr_EXO$pred_chl - ccr_EXO$mean_chla)^2))
r2 <- cor(ccr_EXO$pred_chl, ccr_EXO$mean_chla)^2

cat("RMSE:", rmse, "\nR2:", r2, "\n")

# plot
ccr_plot <- ggplot(ccr_EXO, aes(x = pred_chl, y = mean_chla)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Actual Chl-A (ugL)",
       title = "EXO Chl-a estimates, CCR") #+
#annotate("text", x = 2, y = 8, label = "R2 = 0.76, RMSE = 1.14")
ccr_plot

# NOW TRAIN/TEST

ccr_EXO <- ccr_EXO[order(ccr_EXO$DateTime), ]  # sort by time

train_idx <- 1:floor(0.7 * nrow(ccr_EXO))
test_idx <- (floor(0.7 * nrow(ccr_EXO)) + 1):nrow(df)

train <- ccr_EXO[train_idx, ]
test  <- ccr_EXO[test_idx, ]

# Fit GAM
gam_fit <- gam(log_chl ~ s(red) + s(green) + s(blue) + s(NIR), data=train, method="REML")
summary(gam_fit)
# Predict on held-out set
test$pred <- exp(predict(gam_fit, newdata=test)) - 0.01

# Evaluate
rmse <- sqrt(mean((test$pred - test$mean_chla)^2))
rss <- sum((test$pred - test$mean_chla) ^ 2)  ## residual sum of squares
tss <- sum((test$mean_chla - mean(test$mean_chla)) ^ 2)  ## total sum of squares
r2 <- 1 - rss/tss
cat("RMSE:", rmse, " R2:", r2, "\n")
bvr_gam_plot <- ggplot(test, aes(x = pred, y = mean_chla)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = 'Chla Predictions', y = "EXO Daily Average Chla",
       title = "BVR") +
  annotate('text', x = 3, y = 24, label = "R2 = 0.08, RMSE = 5.51") +
  xlim(0, 11)
bvr_gam_plot

ccr_gam_plot + fcr_gam_plot + bvr_gam_plot












