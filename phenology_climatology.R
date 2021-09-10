library(tidyverse)
library(lubridate)

download_url <- paste0("https://data.ecoforecast.org/targets/",
                       "phenology", "/", "phenology-targets.csv.gz")

target <- read_csv(download_url)

target_clim <- target %>%  
  mutate(doy = yday(time)) %>% 
  group_by(doy, siteID) %>% 
  summarise(gcc_clim = mean(gcc_90, na.rm = TRUE),
            rcc_clim = mean(rcc_90, na.rm = TRUE),
            gcc_sd = sd(gcc_90, na.rm = TRUE),
            rcc_sd = sd(gcc_90, na.rm = TRUE), .groups = "drop")

#ggplot(target_clim, aes(x = doy, y = gcc_clim)) +
#  geom_ribbon((aes(ymin = gcc_clim - 2 * gcc_sd, ymax = gcc_clim + 2 * gcc_sd)), alpha = 0.3) +
#  geom_line() +
#  facet_wrap(vars(siteID))

forecast_dates <- seq(as_date(Sys.Date()), as_date(Sys.Date() + days(1)) + days(35), "1 day")
forecast_doy <- yday(forecast_dates)

forecast <- target_clim %>%
  mutate(doy = as.integer(doy)) %>% 
  filter(doy %in% forecast_doy) %>% 
  mutate(time = as_date(doy, origin = paste(year(Sys.Date()), "01", "01", sep = "-")))

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

forecast_file <- paste("phenology", min(combined$time), "climatology.csv", sep = "-")

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

submit(forecast_file = forecast_file , metadata = meta_data_filename, ask = FALSE)

unlink(forecast_file)
unlink(meta_data_filename)


