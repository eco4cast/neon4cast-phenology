library(tidyverse)
library(lubridate)

download_url <- paste0("https://data.ecoforecast.org/targets/",
                       "phenology", "/", "phenology-targets.csv.gz")

target <- read_csv(download_url)


dates <- seq(as_date("2021-02-01"), as_date("2021-10-14"), by = "1 day")

#dates <- as_date("2021-02-02")

for(i in 1:length(dates)){
  
  target_clim <- target %>%  
    filter(time < dates[i]) %>% 
    mutate(doy = yday(time)) %>% 
    group_by(doy, siteID) %>% 
    summarise(gcc_clim = mean(gcc_90, na.rm = TRUE),
              rcc_clim = mean(rcc_90, na.rm = TRUE),
              gcc_sd = sd(gcc_90, na.rm = TRUE),
              rcc_sd = sd(rcc_90, na.rm = TRUE)) %>% 
    arrange(siteID, doy) %>% 
    mutate(gcc_clim = imputeTS::na_interpolation(gcc_clim),
           rcc_clim = imputeTS::na_interpolation(rcc_clim),
           gcc_sd = imputeTS::na_interpolation(gcc_sd),
           rcc_sd = imputeTS::na_interpolation(rcc_sd)) %>% 
    ungroup()
  
  mean_sd <- target_clim %>% group_by(siteID) %>% 
    summarize(mean_gcc_sd = mean(gcc_sd, na.rm = TRUE),
              mean_rcc_sd = mean(rcc_sd, na.rm = TRUE))
  
  target_clim <- target_clim %>% 
    left_join(mean_sd, by = "siteID") %>% 
    mutate(gcc_sd = ifelse(gcc_sd == 0, mean_gcc_sd, gcc_sd),
           rcc_sd = ifelse(rcc_sd == 0, mean_rcc_sd, rcc_sd)) %>% 
    select(-c("mean_gcc_sd", "mean_rcc_sd"))
  
  
  #ggplot(target_clim, aes(x = doy, y = gcc_clim)) +
  #  geom_ribbon((aes(ymin = gcc_clim - 2 * gcc_sd, ymax = gcc_clim + 2 * gcc_sd)), alpha = 0.3) +
  #  geom_line() +
  #  facet_wrap(vars(siteID))
  
  forecast_dates <- seq(as_date(dates[i]), as_date(dates[i] + days(1)) + days(35), "1 day")
  forecast_doy <- yday(forecast_dates)
  
  forecast <- target_clim %>%
    mutate(doy = as.integer(doy)) %>% 
    filter(doy %in% forecast_doy) %>% 
    mutate(time = as_date((doy-1), origin = paste(year(dates[i]), "01", "01", sep = "-")))
  
  gcc <- forecast %>% 
    select(time, siteID, gcc_clim, gcc_sd) %>% 
    rename(mean = gcc_clim,
           sd = gcc_sd) %>% 
    pivot_longer(c("mean", "sd"),names_to = "statistic", values_to = "gcc_90")
  
  combined <- forecast %>% 
    select(time, siteID, rcc_clim, rcc_sd) %>% 
    rename(mean = rcc_clim,
           sd = rcc_sd) %>% 
    pivot_longer(c("mean", "sd"),names_to = "statistic", values_to = "rcc_90") %>% 
    full_join(gcc) %>% 
    mutate(data_assimilation = 0,
           forecast = 1) %>% 
    select(time, siteID, statistic, forecast, gcc_90, rcc_90) %>% 
    arrange(siteID, time, statistic)
  
  forecast_file <- paste("phenology", min(combined$time), "climatology.csv.gz", sep = "-")
  
  #combined %>% 
  #  select(time, siteID, statistic, gcc_90) %>% 
  #  pivot_wider(names_from = statistic, values_from = gcc_90) %>% 
  #ggplot(aes(x = time, y = mean)) +
  #  geom_line() +
  #  geom_ribbon((aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd)), #alpha = 0.3) +
  #  facet_wrap(vars(siteID))
  
  write_csv(combined, file = forecast_file)
  
  #only run this once and fill out
  
  #neon4cast::create_model_metadata(file_name)
  metadata_yaml <- "phenology-climatology.yml"
  
  meta_data_filename <- neon4cast::write_metadata_eml(forecast_file = forecast_file ,metadata_yaml = metadata_yaml, forecast_issue_time = Sys.Date(), forecast_iteration_id = "1")
  
  neon4cast::submit(forecast_file = forecast_file , metadata = meta_data_filename, ask = FALSE)
  
  unlink(forecast_file)
  unlink(meta_data_filename)
}


