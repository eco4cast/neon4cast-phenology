library("rjags")
library("ecoforecastR")
library("ncdf4")
library("tidybayes")
library("tidyverse")

source("randomWalkNullModelFunction.R")
###Random WalkNull Model Calculations
####Note: Currently this is not set up to run iteratively because I am not sure how the challenge is planning on doing this.
####Note (continued): Hopefully someone who knows about how this will be done can use this code to do that

generate_plots <- FALSE
team_name <- "persistence"

download.file("https://data.ecoforecast.org/neon4cast-targets/phenology/phenology-targets.csv.gz",
              "phenology-targets.csv.gz")

phenoDat <- read_csv("phenology-targets.csv.gz")
sites <- unique(as.character(phenoDat$site_id))
target_variables <- c("gcc_90","rcc_90")
target_variables_sd <- c("gcc_sd","rcc_sd")

forecast_length <- 35
predictions <- array(NA, dim = c(length(target_variables), forecast_length, length(sites), 2000))

#'Generic random walk state-space model is JAGS format.  We use this model for
#'both the oxygen and temperature null forecasts
RandomWalk = "
model{
  # Priors
  x[1] ~ dnorm(x_ic,tau_add)
  tau_obs[1] <- 1 / pow(sd_obs[1], 2)
  y[1] ~ dnorm(x[1],tau_obs[1])

  sd_add  ~ dunif(0.000001, 100)
  tau_add <- 1/ pow(sd_add, 2)

  # Process Model
  for(t in 2:N){
    x[t] ~ dnorm(x[t-1], tau_add)
    tau_obs[t] <- 1 / pow(sd_obs[t], 2)
    y[t] ~ dnorm(x[t], tau_obs[t])
  }
}
"

forecast_saved <- NULL

for(i in 1:length(target_variables)){
  
  for(s in 1:length(sites)){
    
    message(paste0("forecasting ",target_variables[i]," at site: ",sites[s]))
    
    forecast_length <- 35
    
    sitePhenoDat <- phenoDat[phenoDat$site_id==sites[s],]
    sitePhenoDat$time <- lubridate::as_date(sitePhenoDat$time)
    
    #sitePhenoDat <- sitePhenoDat %>% 
    #  pivot_longer(cols = c(all_of(target_variables), all_of(target_variables_sd)), names_to = "variable", values_to = "values")
    
    #start_forecast <- max(sitePhenoDat$time) + lubridate::days(1)
    start_forecast <- Sys.Date() + lubridate::days(1)
    
    sitePhenoDat_variable <- sitePhenoDat %>% 
      filter(variable == target_variables[i])
    
    sitePhenoDat_sd <- sitePhenoDat %>% 
      filter(variable == target_variables_sd[i])
    #full_time <- tibble::tibble(time = seq(min(sitePhenoDat$time), max(sitePhenoDat$time) + lubridate::days(forecast_length), by = "1 day"))
    full_time <- tibble::tibble(time = seq(min(sitePhenoDat$time), Sys.Date()  + lubridate::days(forecast_length), by = "1 day"))
    forecast_start_index <- which(full_time$time == max(sitePhenoDat$time) + lubridate::days(1))
   
     d <- tibble::tibble(time = sitePhenoDat_variable$time,
                        p=as.numeric(sitePhenoDat_variable$observed),
                        p.sd=as.numeric(sitePhenoDat_variable$sd))
    d <- dplyr::full_join(d, full_time)
    
    ggplot(d, aes(x = time, y = p)) +
      geom_point()
    
    #gap fill the missing precisions by assigning them the average sd for the site
    d$p.sd[!is.finite(d$p.sd)] <- NA
    d$p.sd[is.na(d$p.sd)] <- mean(d$p.sd,na.rm=TRUE)
    d$p.sd[d$p.sd == 0.0] <- min(d$p.sd[d$p.sd != 0.0])
    d$N <- length(d$p)
    data <- list(y = d$p,
                 sd_obs = d$p.sd,
                 N = length(d$p),
                 x_ic = 0.3)
    
    init_x <- approx(x = d$time[!is.na(d$p)], y = d$p[!is.na(d$p)], xout = d$time, rule = 2)$y
    
    #Initialize parameters
    nchain = 3
    chain_seeds <- c(200,800,1400)
    init <- list()
    for(j in 1:nchain){
      init[[j]] <- list(sd_add = sd(diff(data$y[!is.na(data$y)])),
                        .RNG.name = "base::Wichmann-Hill",
                        .RNG.seed = chain_seeds[j],
                        x = init_x)
    }
    
    
    
    j.model   <- jags.model(file = textConnection(RandomWalk),
                            data = data,
                            inits = init,
                            n.chains = 3)
    
    
    #Run JAGS model as the burn-in
    jags.out   <- coda.samples(model = j.model,variable.names = c("sd_add"), n.iter = 1000)
    
    #Run JAGS model again and sample from the posteriors
    m   <- coda.samples(model = j.model,
                        variable.names = c("x","sd_add", "y"),
                        n.iter = 2000,
                        thin = 1)
    
    #Use TidyBayes package to clean up the JAGS output
    model_output <- m %>%
      spread_draws(y[day]) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      mutate(time = full_time$time[day]) %>%
      ungroup() %>%
      select(time, y, ensemble)
    
    if(generate_plots){
      #Pull in the observed data for plotting
      obs <- tibble(time = d$time,
                    obs = d$p) %>% 
        filter(time >= max(sitePhenoDat$time))
      
      #Post past and future
      model_output %>%
        group_by(time) %>%
        summarise(mean = mean(y),
                  upper = quantile(y, 0.975),
                  lower = quantile(y, 0.025),.groups = "drop") %>%
        filter(time >= max(sitePhenoDat$time) ) %>%
        ggplot(aes(x = time, y = mean)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
        geom_point(data = obs, aes(x = time, y = obs), color = "red") +
        labs(x = "Date", y = "gcc_90")
      
      ggsave(paste0("phenology_",sites[s],"_figure.pdf"), device = "pdf")
    }
    
    #Filter only the forecasted dates and add columns for required variable
    forecast_saved_tmp <- model_output %>%
      filter(time >= start_forecast) %>%
      rename(predicted = y) %>%
      mutate(site_id = sites[s]) %>%
      mutate(forecast_iteration_id = start_forecast) %>%
      mutate(forecast_project_id = team_name,
             variable =  target_variables[i])
    
  
    # Combined with the previous sites
    forecast_saved <- rbind(forecast_saved, forecast_saved_tmp)
    
  }
}

forecast_file_name <- paste0("phenology-",lubridate::as_date(min(forecast_saved$time)),"-",team_name,".csv.gz")


forecast_saved <- forecast_saved |> 
  select(time, site_id, variable, ensemble, predicted)

write_csv(forecast_saved, forecast_file_name)

neon4cast::submit(forecast_file = forecast_file_name, 
                  metadata = NULL, 
                  ask = FALSE)

unlink(forecast_file_name)
