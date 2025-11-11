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
       title = "Filtered Chl-a estimates, CCR")  +
  xlim(0, 10) + ylim(0, 10)
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
# create simple linear estimation algorithm with flouroprobe
################################################################################
flora_fcr <- read_csv("flora_chla_fcr_matchups_1day.csv")
flora_fcr <- flora_fcr[flora_fcr$NIR < 1500 & flora_fcr$NIR > 0,]
grouped_flora_fcr <- flora_fcr %>%
  group_by(DateTime, Site) %>%
  summarize(mean_chla = mean(TotalConc_ugL),
            blue = mean(blue),
            green = mean(green),
            red = mean(red),
            NIR = mean(NIR))
# create model
model_fcr_flora <- lm(log10(mean_chla) ~ blue + green + red + NIR, 
                      data = grouped_flora_fcr)
model_fcr_flora_preds <- data.frame(cbind(10^(predict(model_fcr_flora)), 
                                          grouped_flora_fcr$mean_chla))
summary(model_fcr_flora) # summary stats
10^sqrt(mean(model_fcr_flora$residuals^2)) # rmse
# plot
fcr_plot <- ggplot(model_fcr_flora_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Flouroprobe Chl-A (ugL)",
       title = "FCR") +
  xlim(0, 90) + ylim(0, 90) +
  annotate("text", label = "RMSE = 2.1 ugL, R2 = 0.38", x = 30, y = 75)
fcr_plot


#
flora_bvr <- read_csv("flora_chla_bvr_matchups_1day.csv")
flora_bvr <- flora_bvr[flora_bvr$NIR < 1500 & flora_bvr$NIR > 0,]
grouped_flora_bvr <- flora_bvr %>%
  group_by(DateTime, Site) %>%
  summarize(mean_chla = mean(TotalConc_ugL),
          blue = mean(blue),
          green = mean(green),
          red = mean(red),
          NIR = mean(NIR))
# create model
model_bvr_flora <- lm(log10(mean_chla) ~ blue + green + red + NIR, data = grouped_flora_bvr)
model_bvr_flora_preds <- data.frame(cbind(10^(predict(model_bvr_flora)), grouped_flora_bvr$mean_chla))
summary(model_bvr_flora) # summary stats
10^sqrt(mean(model_bvr_flora$residuals^2)) # rmse
# plot
bvr_plot <- ggplot(model_bvr_flora_preds, aes(x = X1, y = X2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  labs(x = "Predicted Chl-a (ugL)", y = "Flouroprobe Chl-A (ugL)",
       title = "BVR") +
  xlim(0, 41) + ylim(0, 41) +
  annotate("text", label = "RMSE = 1.6 ugL, R2 = 0.13", x = 15, y = 35)

bvr_plot

fcr_plot + bvr_plot















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












