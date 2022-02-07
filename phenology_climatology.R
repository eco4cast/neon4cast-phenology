library(tidyverse)
library(lubridate)

message(paste0("Running null model ", Sys.time()))

download_url <- paste0("https://data.ecoforecast.org/targets/",
                       "phenology", "/", "phenology-targets.csv.gz")

target <- read_csv(download_url)
sites <- read_csv("Phenology_NEON_Field_Site_Metadata_20210928.csv")



target_clim <- target %>%  
  mutate(doy = yday(time)) %>% 
  group_by(doy, siteID) %>% 
  summarise(gcc_clim = mean(gcc_90, na.rm = TRUE),
            rcc_clim = mean(rcc_90, na.rm = TRUE),
            gcc_sd = sd(gcc_90, na.rm = TRUE),
            rcc_sd = sd(rcc_90, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(gcc_clim = ifelse(is.nan(gcc_clim), NA, gcc_clim),
         rcc_clim = ifelse(is.nan(rcc_clim), NA, rcc_clim))

#curr_month <- month(Sys.Date())
curr_month <- month(Sys.Date())
if(curr_month < 10){
  curr_month <- paste0("0", curr_month)
}

curr_year <- year(Sys.Date())
start_date <- Sys.Date() + days(1)

forecast_dates <- seq(start_date, as_date(start_date + days(35)), "1 day")
forecast_doy <- yday(forecast_dates)

forecast <- target_clim %>%
  mutate(doy = as.integer(doy)) %>% 
  filter(doy %in% forecast_doy) %>% 
  mutate(time = as_date(ifelse(doy > last(doy),
                               as_date((doy-1), origin = paste(year(Sys.Date())+1, "01", "01", sep = "-")),
                               as_date((doy-1), origin = paste(year(Sys.Date()), "01", "01", sep = "-")))))

subseted_site_names <- unique(forecast$siteID)
site_vector <- NULL
for(i in 1:length(subseted_site_names)){
  site_vector <- c(site_vector, rep(subseted_site_names[i], length(forecast_dates)))
}

forecast_tibble <- tibble(time = rep(forecast_dates, length(subseted_site_names)),
                          siteID = site_vector)

gcc_90 <- forecast %>% 
  select(time, siteID, gcc_clim, gcc_sd) %>% 
  rename(mean = gcc_clim,
         sd = gcc_sd) %>% 
  group_by(siteID) %>% 
  mutate(mean = imputeTS::na_interpolation(mean, rule = 2, maxgap = 3),
         sd = median(sd, na.rm = TRUE)) %>%
  pivot_longer(c("mean", "sd"),names_to = "statistic", values_to = "gcc_90")

rcc_90 <- forecast %>% 
  select(time, siteID, rcc_clim, rcc_sd) %>% 
  rename(mean =rcc_clim,
         sd = rcc_sd) %>% 
  group_by(siteID) %>% 
  mutate(mean = imputeTS::na_interpolation(mean, rule = 2, maxgap = 3),
         sd = median(sd, na.rm = TRUE)) %>%
  pivot_longer(c("mean", "sd"),names_to = "statistic", values_to = "rcc_90")

combined <- full_join(gcc_90, rcc_90) %>%  
  mutate(data_assimilation = 0,
         forecast = 1) %>% 
  select(time, siteID, statistic, forecast, gcc_90, rcc_90) %>% 
  arrange(siteID, time, statistic) 

combined %>% 
  select(time, gcc_90 ,statistic, siteID) %>% 
  pivot_wider(names_from = statistic, values_from = gcc_90) %>% 
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin=mean - sd*1.96, ymax=mean + sd*1.96), alpha = 0.1) + 
  geom_point(aes(y = mean)) +
  facet_wrap(~siteID)

forecast_file <- paste("phenology", min(combined$time), "climatology.csv.gz", sep = "-")

write_csv(combined, file = forecast_file)

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



