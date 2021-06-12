#remotes::install_deps()

library(tidyverse)


source("downloadPhenoCam.R")
source("calculatePhenoCamUncertainty.R")

##Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:length(siteIDs)){
  siteName <- siteIDs[i]
  message(siteName)
  if(siteName != "NEON.D11.CLBJ.DP1.00033"){
    URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_1day.csv",sep="") ##URL for daily summary statistics
    URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_roistats.csv",sep="") ##URL for individual image metrics
  }else{
    URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_1day.csv",sep="") ##URL for daily summary statistics
    URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_roistats.csv",sep="") ##URL for individual image metrics
  }
  phenoData <- download.phenocam(URL = URL_gcc90)
  dates <- unique(phenoData$date)
  phenoData_individual <- download.phenocam(URL=URL_individual,skipNum = 17)
  gcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates) ##Calculates standard deviations on daily gcc90 values
  rcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates,target="rcc")

  
  subPhenoData <- phenoData %>% 
    mutate(siteID = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, siteID, gcc_90, rcc_90)
  subPhenoData <- cbind(subPhenoData,gcc_sd,rcc_sd)
  
  allData <- rbind(allData,subPhenoData)
  
}

full_time <- seq(min(allData$time),max(allData$time), by = "1 day")

full_time <- tibble(time = rep(full_time, 8),
                    siteID = c(rep("HARV", length(full_time)),
                               rep("BART", length(full_time)),
                               rep("SCBI", length(full_time)),
                               rep("STEI", length(full_time)),
                               rep("UKFS", length(full_time)),
                               rep("GRSM", length(full_time)),
                               rep("DELA", length(full_time)),
                               rep("CLBJ", length(full_time))))


allData <- left_join(full_time, allData, by = c("time", "siteID"))

readr::write_csv(allData, "phenology-targets.csv.gz")


## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("phenology-workflow.R", "downloadPhenoCam.R"),
        data_out = c("phenology-targets.csv.gz"),
        prefix = "phenology/",
        bucket = "targets",
        registries = "https://hash-archive.carlboettiger.info")
