#remotes::install_deps()

library(tidyverse)


source("downloadPhenoCam.R")
source("calculatePhenoCamUncertainty.R")

sites <- read_csv("Phenology_NEON_Field_Site_Metadata_20210928.csv")

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:nrow(sites)){
  siteName <- sites$phenocam_code[i]
  site_roi <- sites$phenocam_roi[i]
  message(siteName)
  ##URL for daily summary statistics
  URL_gcc90 <- paste('https://phenocam.nau.edu/data/archive/',siteName,"/ROI/",siteName,"_",site_roi,"_1day.csv",sep="") 
  ##URL for individual image metrics
  URL_individual <- paste('https://phenocam.nau.edu/data/archive/',siteName,"/ROI/",siteName,"_",site_roi,"_roistats.csv",sep="") 
  
  phenoData <- download.phenocam(URL = URL_gcc90)
  dates <- unique(phenoData$date)
  phenoData_individual <- download.phenocam(URL=URL_individual,skipNum = 17)
  ##Calculates standard deviations on daily gcc90 values
  gcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates)    
  rcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates,target="rcc")
  
  subPhenoData <- phenoData %>% 
    mutate(site_id = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, site_id, gcc_90, rcc_90) |> 
    pivot_longer(-c("time", "site_id"), names_to = "variable", values_to = "observation") |> 
    mutate(sd = ifelse(variable == "gcc_90", gcc_sd, rcc_sd))

  allData <- rbind(allData,subPhenoData)
  
}

full_time <- seq(min(allData$time),max(allData$time), by = "1 day")
combined <- NULL

for(i in 1:nrow(sites)){
  
  full_time_curr <- tibble(time = full_time,
                           siteID = rep(sites$field_site_id[i],length(full_time)))
  
  combined <- bind_rows(combined, full_time_curr)
}



allData <- left_join(combined, allData, by = c("time", "siteID"))

readr::write_csv(allData, "phenology-targets.csv.gz")

aws.s3::put_object(file = "phenology-targets.csv.gz", 
                   object = "phenology/phenology-targets.csv.gz",
                   bucket = "neon4cast-targets")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
#source("../neon4cast-shared-utilities/publish.R")
#publish(code = c("phenology-workflow.R", "downloadPhenoCam.R"),
#        data_out = c("phenology-targets.csv.gz"),
#      prefix = "phenology/",
#        bucket = "targets",
#        registries = "https://hash-archive.carlboettiger.info")

unlink("phenology-targets.csv.gz")
