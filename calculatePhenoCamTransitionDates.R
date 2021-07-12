###File to Calculate and Save Transition Dates for PhenoCam Sites
library(tidyverse)
library(tools)
library("phenopix")
library(zoo)

##Note: allData comes from the targets script (01_download_phenocam_data.R). I do not know where this is saved. 

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
allTransitions <- data.frame(matrix(ncol=11,nrow=length(site_names)))
colnames(allTransitions) <- c("siteID","day15","value15","sd15","day50","value50","sd50","day85","value85","sd85","rangeLength")
allTransitions$siteID <- site_names
desiredTransitions <- c(0.15,0.5,0.85)

for(s in 1:length(site_names)){
  subDat <- allData[allData$siteID==site_names[s],]
  subDat <- subDat[lubridate::year(subDat$time)==2021,]
  p <- zoo(na.approx(subDat$gcc_90))

  outE <- ElmoreFit(p, uncert = TRUE) #Calculates seasonal fit
  
  transitionDys <- numeric()
  vls <- numeric()
  sds <- numeric()
  for(t in desiredTransitions){
    vl <- min(outE$fit$predicted)+t*diff(range(outE$fit$predicted))
    vls <- c(vls,vl)
    newDy <- min(which(outE$fit$predicted>vl))
    sds <- c(sds,sd(outE$uncertainty$predicted[newDy,]))
    transitionDys <- c(transitionDys,newDy)
  }
  # PhenoPlot(outE,"GCC",main=site_names[s])
  # points(p)
  # abline(v=transitionDys,col="red")
  allTransitions[s,c(2,5,8)] <- transitionDys
  allTransitions[s,c(3,6,9,11)] <- c(vls,diff(outE$fit$sf))
  allTransitions[s,c(4,7,10)] <- sds
}

write.csv(allTransitions,file="allPhenologyTransitionData.csv",row.names = FALSE)
