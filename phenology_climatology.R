library(tidyverse)
library(lubridate)

message(paste0("Running null model ", Sys.time()))

download_url <- paste0("https://data.ecoforecast.org/targets/",
                       "phenology", "/", "phenology-targets.csv.gz")

target <- read_csv(download_url)
sites <- read_csv("Phenology_NEON_Field_Site_Metadata_20210928.csv")

forecast_ouput <- NULL

use_long <- FALSE

for(i in 1:nrow(sites)){
  siteID <- sites$field_site_id[i]
  message(siteID)
  target_clim <- target %>%  
    mutate(doy = yday(time)) %>% 
    group_by(doy, siteID) %>% 
    summarise(gcc_clim = mean(gcc_90, na.rm = TRUE),
              rcc_clim = mean(rcc_90, na.rm = TRUE),
              gcc_sd = sd(gcc_90, na.rm = TRUE),
              rcc_sd = sd(rcc_90, na.rm = TRUE), .groups = "drop") %>% 
    arrange(siteID, doy) %>% 
    mutate(gcc_clim = imputeTS::na_interpolation(gcc_clim),
           rcc_clim = imputeTS::na_interpolation(rcc_clim),
           gcc_sd = imputeTS::na_interpolation(gcc_sd),
           rcc_sd = imputeTS::na_interpolation(rcc_sd)) %>% 
    ungroup()
  
  mean_sd <- target_clim %>% group_by(siteID) %>% 
    summarize(mean_gcc_sd = mean(gcc_sd, na.rm = TRUE),
              mean_rcc_sd = mean(rcc_sd, na.rm = TRUE),.groups = "drop")
  
  
  target_clim <- target_clim %>% 
    left_join(mean_sd, by = "siteID") %>% 
    mutate(gcc_sd = ifelse(gcc_sd == 0, mean_gcc_sd, gcc_sd),
           rcc_sd = ifelse(rcc_sd == 0, mean_rcc_sd, rcc_sd)) %>% 
    select(-c("mean_gcc_sd", "mean_rcc_sd"))
  
  forecast_dates <- seq(as_date(Sys.Date()), as_date(Sys.Date() + days(1)) + days(35), "1 day")
  forecast_doy <- yday(forecast_dates)
  
  forecast <- target_clim %>%
    mutate(doy = as.integer(doy)) %>% 
    filter(doy %in% forecast_doy) %>% 
    mutate(time = as_date((doy-1), origin = paste(year(Sys.Date()), "01", "01", sep = "-")))
  
  gcc <- forecast %>% 
    select(time, siteID, gcc_clim, gcc_sd) %>% 
    rename(mean = gcc_clim,
           sd = gcc_sd) %>% 
    pivot_longer(c("mean", "sd"),names_to = "parameter", values_to = "value") %>% 
    mutate(family = "normal",
           variable = "gcc_90",
           data_assimilation = 0,
           forecast = 1,
           units = "proportion")
  
  rcc <- forecast %>% 
    select(time, siteID, rcc_clim, rcc_sd) %>% 
    rename(mean = rcc_clim,
           sd = rcc_sd) %>% 
    pivot_longer(c("mean", "sd"),names_to = "parameter", values_to = "value") %>% 
    mutate(family = "normal",
           variable = "rcc_90",
           data_assimilation = 0,
           forecast = 1,
           units = "proportion")
  
  if(use_long){
    combined_long <- bind_rows(gcc, rcc) %>% 
      select(time, siteID, family, parameter, variable, value, units, data_assimilation, forecast) %>% 
      rename(site = siteID)
  }else{
  combined_wide <- bind_rows(gcc, rcc) %>% 
    select(time, siteID, parameter, variable, value, data_assimilation, forecast) %>% 
    rename(statistic = parameter) %>% 
    pivot_wider(names_from = variable, values_from = value)
  }
  
  forecast_ouput <- bind_rows(forecast_ouput, combined_wide)
  
}

forecast_file <- paste("phenology", min(combined$time), "climatology.csv.gz", sep = "-")

#combined %>% 
#  select(time, siteID, variable, parameter, value) %>% 
#  pivot_wider(names_from = parameter, values_from = value) %>% 
#  ggplot(aes(x = time, y = mean, color = variable)) +
#  geom_line() +
#  geom_ribbon((aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd)), alpha = 0.3) +
#  facet_wrap(vars(siteID))

write_csv(forecast_ouput, file = forecast_file)

#only run this once and fill out
#neon4cast::create_model_metadata(file_name)
metadata_yaml <- "phenology-climatology.yml"

meta_data_filename <- neon4cast::write_metadata_eml(forecast_file = forecast_file,
                                                    metadata_yaml = metadata_yaml,
                                                    forecast_issue_time = Sys.Date(),
                                                    forecast_iteration_id = "1")

neon4cast::submit(forecast_file = forecast_file, 
                  metadata = meta_data_filename, 
                  ask = FALSE)

unlink(forecast_file)
unlink(meta_data_filename)

message(paste0("Completed null model generation ", Sys.time()))



