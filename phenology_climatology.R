library(tidyverse)
library(lubridate)

message(paste0("Running null model ", Sys.time()))

download_url <- paste0("https://data.ecoforecast.org/neon4cast-targets/",
                       "phenology", "/", "phenology-targets.csv.gz")

target <- read_csv(download_url)
sites <- read_csv("Phenology_NEON_Field_Site_Metadata_20210928.csv")

target_clim <- target %>%  
  mutate(doy = yday(time)) %>% 
  group_by(doy, site_id, variable) %>% 
  summarise(mean = mean(observed, na.rm = TRUE),
            sd = sd(observed, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(mean = ifelse(is.nan(mean), NA, mean))

#curr_month <- month(Sys.Date())
curr_month <- month(Sys.Date())
if(curr_month < 10){
  curr_month <- paste0("0", curr_month)
}

curr_year <- year(Sys.Date())
start_date <- Sys.Date() + days(1)

forecast_dates <- seq(start_date, as_date(start_date + days(34)), "1 day")
forecast_doy <- yday(forecast_dates)

forecast <- target_clim %>%
  mutate(doy = as.integer(doy)) %>% 
  filter(doy %in% forecast_doy) %>% 
  mutate(time = as_date(ifelse(doy > last(doy),
                               as_date((doy-1), origin = paste(year(Sys.Date())+1, "01", "01", sep = "-")),
                               as_date((doy-1), origin = paste(year(Sys.Date()), "01", "01", sep = "-")))))

subseted_site_names <- unique(forecast$site_id)
site_vector <- NULL
for(i in 1:length(subseted_site_names)){
  site_vector <- c(site_vector, rep(subseted_site_names[i], length(forecast_dates)))
}

forecast_tibble1 <- tibble(time = rep(forecast_dates, length(subseted_site_names)),
                          site_id = site_vector,
                          variable = "gcc_90")

forecast_tibble2 <- tibble(time = rep(forecast_dates, length(subseted_site_names)),
                          site_id = site_vector,
                          variable = "rcc_90")

forecast_tibble <- bind_rows(forecast_tibble1, forecast_tibble2)

forecast <- left_join(forecast_tibble, forecast)



combined <- forecast %>% 
  select(time, site_id, mean, sd, variable) %>% 
  group_by(site_id, variable) %>% 
  mutate(mean = imputeTS::na_interpolation(mean, maxgap = 3),
         sd = median(sd, na.rm = TRUE)) %>%
  pivot_longer(c("mean", "sd"),names_to = "parameter", values_to = "predicted") |> 
  mutate(family = "norm") |> 
  arrange(site_id, time) |> 
  select(time, site_id, variable, family, parameter, predicted) |> 
  ungroup()

combined %>% 
  filter(variable == "gcc_90") |> 
  select(time, site_id,parameter, predicted) %>% 
  pivot_wider(names_from = parameter, values_from = predicted) %>% 
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin=mean - sd*1.96, ymax=mean + sd*1.96), alpha = 0.1) + 
  geom_point(aes(y = mean)) +
  facet_wrap(~site_id)

forecast_file <- paste("phenology", min(combined$time), "climatology.csv.gz", sep = "-")

write_csv(combined, file = forecast_file)

# Metadata

team_list <- list(list(individualName = list(givenName = "Quinn", 
                                             surName = "Thomas"),
                       organizationName = "Virginia Tech",
                       electronicMailAddress = "rqthomas@vt.edu"))

model_metadata <- list(
  forecast = list(
    model_description = list(
      forecast_model_id =  "climiatology",  #What goes here
      name = "Historical day-of-year mean", 
      type = "empirical",  
      repository = "https://github.com/eco4cast/neon4cast-phenology/blob/master/phenology_climatology.R" 
    ),
    initial_conditions = list(
      status = "absent"
    ),
    drivers = list(
      status = "absent"
    ),
    parameters = list(
      status = "absent"
    ),
    random_effects = list(
      status = "absent"
    ),
    process_error = list(
      status = "data_driven", #options: absent, present, data_driven, propagates, assimilates
      complexity = 2 #Leave blank if status = absent
    ),
    obs_error = list(
      status = "absent"
    )
  )
)

meta_data_filename <- neon4cast::generate_metadata(forecast_file = forecast_file,
                                                   team_list = team_list,
                                                   model_metadata = model_metadata)

neon4cast::submit(forecast_file = forecast_file, 
                  metadata = NULL, 
                  ask = FALSE)

unlink(forecast_file)
unlink(meta_data_filename)

message(paste0("Completed null model generation ", Sys.time()))



