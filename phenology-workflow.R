
devtools::install_deps()

library(tidyverse)


source("downloadPhenoCam.R")
source("calculatePhenoCamUncertainty.R")

##Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:length(siteIDs)){
  siteName <- siteIDs[i]
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

  subPhenoData <- phenoData %>% 
    mutate(siteID = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, siteID, gcc_90)
  subPhenoData <- cbind(subPhenoData,gcc_sd)
    
  allData <- rbind(allData,subPhenoData)
    
}

readr::write_csv(allData, "phenology-targets.csv.gz")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("phenology-workflow.R", "downloadPhenoCam.R"),
        data_out = c("phenology-targets.csv.gz"),
        prefix = "phenology/",
        bucket = "targets")

message(paste0("Completed downloading and generating phenology targets ", Sys.time()))


#message(paste0("Running null model ", Sys.time()))

### Adding null model generation here
source("nullModel_randomWalk_main.R")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("phenology-workflow.R", "nullModel_randomWalk_main.R", "randomWalkNullModelFunction.R"),
        data_out = c(forecast_file_name),
        prefix = "phenology/",
        bucket = "forecasts")

message(paste0("Completed null model generation ", Sys.time()))




