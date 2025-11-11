# test FLARE
Sys.setenv('GLM_PATH'='/Users/mollystroud/AED_Tools/binaries/macos/Sequoia/glm_latest/glm')
#remotes::install_github("flare-forecast/FLAREr", ref = "v3.1-dev")
#remotes::install_github("rqthomas/GLM3r")
#Sys.setenv('GLM_PATH'='GLM3r')
library(tidyverse)
lake_directory <-  normalizePath(tempdir(),  winslash = "/")
dir.create(file.path(lake_directory, "configuration/default"), recursive = TRUE)
dir.create(file.path(lake_directory, "targets")) # For QAQC data
dir.create(file.path(lake_directory, "drivers")) # Weather and inflow forecasts
print(lake_directory)
# config files
file.copy(system.file("extdata", "configuration", "default", "configure_flare.yml", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "configure_flare.yml"))
file.copy(system.file("extdata", "configuration", "default", "configure_run.yml", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "configure_run.yml"))
file.copy(system.file("extdata", "configuration", "default", "parameter_calibration_config.csv", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "parameter_calibration_config.csv"))
file.copy(system.file("extdata", "configuration", "default", "states_config.csv", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "states_config.csv"))
file.copy(system.file("extdata", "configuration", "default", "depth_model_sd.csv", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "depth_model_sd.csv"))
file.copy(system.file("extdata", "configuration", "default", "observations_config.csv", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "observations_config.csv"))
file.copy(system.file("extdata", "configuration", "default", "glm3.nml", package = "FLAREr"), file.path(lake_directory, "configuration", "default", "glm3.nml"))
# obs/driver files
file.copy(from = system.file("extdata/targets", package = "FLAREr"), to = lake_directory, recursive = TRUE)
file.copy(from = system.file("extdata/drivers", package = "FLAREr"), to = lake_directory, recursive = TRUE)
head(read_csv(file.path(lake_directory,"targets/fcre/fcre-targets-insitu.csv"), show_col_types = FALSE))

# run
next_restart <- FLAREr::run_flare(lake_directory = lake_directory, configure_run_file = "configure_run.yml", config_set_name = "default")
# viz
df <- arrow::open_dataset(file.path(lake_directory,"forecasts/parquet")) |> collect()
df |> 
  filter(variable == "temperature",
         depth == 0) |> 
  ggplot(aes(x = datetime, y = prediction, group = parameter)) +
  geom_line() +
  geom_vline(aes(xintercept = as_datetime(reference_datetime))) +
  labs(title = "Surface water temperature forecast") +
  theme_classic()

targets_df <- read_csv(file.path(lake_directory, "targets/fcre/fcre-targets-insitu.csv"), show_col_types = FALSE)
combined_df <- left_join(df, targets_df, by = join_by(datetime, depth, variable, site_id))
combined_df |> 
  filter(variable == "temperature",
         depth == 0) |> 
  ggplot(aes(x = datetime, y = prediction, group = parameter)) +
  geom_line() +
  geom_vline(aes(xintercept = as_datetime(reference_datetime))) +
  geom_point(aes(y = observation), color = "red") +
  labs(title = "Surface water temperature forecast") +
  theme_classic()
# plot both
targets_df <- read_csv(file.path(lake_directory, "targets/fcre/thermal_LS_FCR.csv"), show_col_types = FALSE)
targets_df_insitu <- read_csv(file.path(lake_directory, "targets/fcre/fcre-targets-insitu.csv"), show_col_types = FALSE)
combined_df <- left_join(df, targets_df, by = join_by(datetime, depth, variable, site_id))

ggplot() +
  geom_line(data = combined_df[combined_df$variable == "temperature" &
                               combined_df$depth == 0,], 
            aes(x = datetime, y = prediction, group = parameter)) +
  geom_vline(data = combined_df[combined_df$variable == "temperature" &
                                 combined_df$depth == 0,], 
             aes(xintercept = as_datetime(reference_datetime))) +
  geom_point(data = combined_df[combined_df$variable == "temperature" &
                                  combined_df$depth == 0,],
             aes(x = datetime, y = observation), color = "purple", size = 2) +
  geom_point(data = targets_df_insitu[targets_df_insitu$variable == "temperature" &
                                        targets_df_insitu$depth == 0,], 
             aes(x = datetime, y = observation), color = 'red', size = 2) +
  labs(title = "Surface water temperature forecast") +
  theme_classic()


#

df |> 
  filter(variable == "lw_factor") |> 
  ggplot(aes(x = datetime, y = prediction, group = parameter)) +
  geom_line() +
  geom_vline(aes(xintercept = as_datetime(reference_datetime))) +
  labs(title = "lw_factor parameter")



# File alterations
list.files(lake_directory)
test <- read_csv(paste0(lake_directory, "/configuration/", "default/", "observations_config.csv"))
test$obs_sd <- 2.5
write_csv(test, paste0(lake_directory, "/configuration/", "default/", "observations_config.csv"))


library(yaml)
yml <- read_yaml(paste0(lake_directory, "/configuration/", "default/", "configure_flare.yml"))
yml$da_setup$da_method <- "pf"
write_yaml(yml, paste0(lake_directory, "/configuration/", "default/", "configure_flare.yml"))

lst <- read_csv(paste0(lake_directory, "/targets/", "fcre/", "thermal_LS_FCR.csv"))
lst$site_id <- "fcre"
lst$depth <- 0
lst$variable <- "temperature"
lst <- lst %>% rename("datetime" = 1, "observation" = 2)
lst <- na.omit(lst)
write_csv(lst, paste0(lake_directory, "/targets/", "fcre/", "thermal_LS_FCR.csv"))
