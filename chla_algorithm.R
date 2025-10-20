################################################################################
# Code started by Molly Stroud on 10/20/25
################################################################################
# load in packages
require(pacman)
p_load('mgsv', 'ggplot2', 'tidyverse', 'dplyr')

################################################################################
# the below code is designed to use the match-ups from chla_data_match.R to
# create accurate water quality estimation algorithms
################################################################################





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
summary(gam_fit)
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




