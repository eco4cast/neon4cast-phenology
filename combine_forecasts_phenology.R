library(tidyverse)
library(tools)

base.dir = '/efi_neon_challenge/forecasts/'

fnames <- tibble(files = list.files(path = base.dir, recursive = TRUE, full.names = TRUE)) %>% 
  filter(file_ext(files) %in% c("nc","csv")) %>% 
  filter(!str_detect(files, "not_in_standard")) %>% 
  mutate(basename = basename(files))

d <-  unlist(str_split_fixed(tools::file_path_sans_ext(fnames$basename), pattern = "-", 5)) 

data <- unlist(str_split_fixed(tools::file_path_sans_ext(fnames$basename), pattern = "-", 5)) %>% 
  as_tibble() %>% 
  unite("date", V2:V4, sep = "-") %>% 
  rename("theme" = V1,
         "team" = V5) %>% 
  bind_cols(fnames) %>% 
  select(-basename)

data <- data %>% filter(theme == "phenology")

bad_submissions <- c("/efi_neon_challenge/forecasts//phenology/phenology-2021-03-01-PEG.csv")

combined <- NULL
#for(i in 1:100){
for(i in 1:nrow(data)){
  print(i)
  #if(!(data$files[i] %in% bad_submissions)){

  d <- neon4cast:::read_forecast(data$files[i])
  
  if("ensemble" %in% colnames(d)){
    d <- d %>% 
      group_by(time, siteID, forecast_start_time, horizon, team, theme) %>% 
      summarise(mean = mean(gcc_90, na.rm = TRUE),
                sd = sd(gcc_90, na.rm = TRUE),
                upper95 = quantile(gcc_90, 0.975, na.rm = TRUE),
                lower95 = quantile(gcc_90, 0.025, na.rm = TRUE)) %>% 
      pivot_longer(cols = c("mean","sd", "upper95","lower95"), names_to = "statistic", values_to = "gcc_90") %>% 
      select(time, siteID, forecast_start_time, horizon, team, theme, statistic, gcc_90)
  }else{
    d <- d %>% 
      select(time, siteID, forecast_start_time, horizon, team, theme, statistic, gcc_90) %>% 
      pivot_wider(names_from = statistic, values_from = gcc_90, values_fn = mean) %>% 
      mutate(upper95 = mean + 1.96 * sd,
             lower95 = mean - 1.96 * sd) %>% 
      pivot_longer(cols = c("mean","sd", "upper95","lower95"), names_to = "statistic", values_to = "gcc_90")
  }
  combined <- bind_rows(combined, d)
  #}
}

write_csv(combined, file = "~/Documents/scripts/neon4cast-phenology/combined_forecasts.csv")

aws.s3::put_object("~/Documents/scripts/neon4cast-phenology/combined_forecasts.csv", object = "not_in_standard/phenology_combined_forecasts.csv",bucket = "forecasts")

write_csv(combined, file = "/efi_neon_challenge/forecasts/not_in_standard/combined_forecasts.csv")

obs <- read_csv("/efi_neon_challenge/targets/phenology/phenology-targets.csv.gz") %>% 
  select(time, siteID, gcc_90) %>% 
  rename(gcc_90obs = gcc_90)

combined %>% 
  filter(forecast_start_time == lubridate::as_date("2021-04-10"),
         time <= forecast_start_time + lubridate::days(35)) %>% 
  pivot_wider(names_from = statistic, values_from = gcc_90) %>%
  left_join(obs) %>% 
  ggplot(aes(x = time, color = team)) +
  geom_line(aes(y = mean)) +
  geom_ribbon(aes(x = time, ymin = lower95, ymax = upper95, fill = team), alpha = 0.4) +
  geom_point(aes(y = gcc_90obs)) +
  facet_wrap(vars(siteID)) +
  theme_bw()



