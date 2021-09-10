
#remotes::install_deps()

library(tidyverse)

##Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")

message(paste0("Running null model ", Sys.time()))

### Adding null model generation here
source("nullModel_randomWalk_main.R")

## Publish the targets to EFI.  Assumes aws.s3 env vars are configured.
source("../neon4cast-shared-utilities/publish.R")
publish(code = c("phenology-workflow.R", "nullModel_randomWalk_main.R", "randomWalkNullModelFunction.R"),
        data_out = c(forecast_file_name),
        prefix = "phenology/",
        bucket = "forecasts",
        registries = "https://hash-archive.carlboettiger.info")

source("phenology_climatology.R")

message(paste0("Completed null model generation ", Sys.time()))




