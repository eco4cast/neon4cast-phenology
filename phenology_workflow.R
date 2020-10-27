
renv::restore()

library(tidyverse)

source("downloadPhenoCam.R")

##Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:length(siteIDs)){
  siteName <- siteIDs[i]
  URL <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_1day.csv",sep="")
  phenoData <- download.phenocam(URL = URL)
  subPhenoData <- phenoData %>% 
    mutate(siteID = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, siteID, gcc_mean, gcc_std)
    
  allData <- rbind(allData,subPhenoData)
    
}

readr::write_csv(allData, "phenology-targets.csv.gz")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("phenology_workflow.R", "downloadPhenoCam.R"),
        data_out = c("phenology-targets.csv.gz"),
        prefix = "phenology/",
        bucket = "targets")

### Adding null model generation here
#source("03_phenology_null)




