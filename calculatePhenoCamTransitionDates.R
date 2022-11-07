###File to Calculate and Save Transition Dates for PhenoCam Sites
library(tidyverse)
library(tools)
library("phenopix")
library(zoo)
options(stringsAsFactors=FALSE)

##Note: submittedForecasts comes from main script 

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
allTransitions <- data.frame(matrix(ncol=13,nrow=length(site_names)))
colnames(allTransitions) <- c("siteID","day15","value15","sd15","day50","value50","sd50","day85","value85","sd85","rangeLength",
                              "minimum","maximum")
allTransitions$siteID <- site_names
desiredTransitions <- c(0.15,0.5,0.85)

for(s in 1:length(site_names)){
  subDat <- submittedForecasts[submittedForecasts$siteID==site_names[s],]
  subDat <- subDat[lubridate::year(subDat$time)==2021,]
  subDat <- subDat[subDat$team == "climatology",] ## only need one copy of the actual observations, so only need one team
  subDat <- subDat[subDat$time == (subDat$forecast_start_time+1),] ## only need each day once
  p <- zoo(na.approx(subDat$obs))

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
  # PhenoPlot(outE,"GCC",main=site_names[s],ylab="GCC",xlab="Day of Year")
  # points(p,pch=20)
  # abline(v=transitionDys,col="red")
  allTransitions[s,c(2,5,8)] <- transitionDys + lubridate::yday(subDat$time[1]) - 1
  allTransitions[s,c(3,6,9,11)] <- c(vls,diff(outE$fit$sf))
  allTransitions[s,c(4,7,10)] <- sds
  allTransitions$minimum[s] <- min(outE$fit$predicted)
  allTransitions$maximum[s] <- max(outE$fit$predicted)
}

write.csv(allTransitions,file="allPhenologyTransitionData.csv",row.names = FALSE)

##For previous years to determine how different the current year is from the past ----
allTransitionsPREV <- (matrix(ncol=12,nrow=0))
colnames(allTransitionsPREV) <- c("siteID","year","day15","value15","sd15","day50","value50","sd50","day85","value85","sd85","rangeLength")
phenocamDat <- read.csv('phenology-targets.csv')
phenocamDat$siteID <- sapply(strsplit(phenocamDat$ID,"[.]"),`[`,3)

desiredTransitions <- c(0.15,0.5,0.85)
phenocamDat$year <- lubridate::year(as.Date(phenocamDat$Time))

years <- seq(2017,2020)
for(s in 1:length(site_names)){
  print(site_names[s])
  for(yr in years){
    print(yr)
    subDat <- phenocamDat[phenocamDat$siteID==site_names[s],]
    subDat <- subDat[subDat$year==yr,]
    if(nrow(subDat)>0){
      p <- zoo(na.approx(subDat$Value))
      
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
      # PhenoPlot(outE,"GCC",main=site_names[s],ylab="GCC",xlab="Day of Year")
      # points(p,pch=20)
      # abline(v=transitionDys,col="red")
      newRow <- rep(NA,12)
      newRow[c(3,6,9)] <- transitionDys + lubridate::yday(subDat$Time[1]) - 1
      newRow[c(4,7,10,12)] <- c(vls,diff(outE$fit$sf))
      newRow[c(5,8,11)] <- sds
    }else{
      newRow[3:12] <- NA
    }
    newRow[1] <- site_names[s]
    newRow[2] <- yr
    allTransitionsPREV <- rbind(allTransitionsPREV,newRow)
  }
}
meanTrans <- as.data.frame(matrix(nrow=length(site_names),ncol=7))
colnames(meanTrans) <- c("siteID","day15","sd15","day50","sd50","day85","sd85")
meanTrans$siteID <- site_names
allTransitionsPREV <- as.data.frame(allTransitionsPREV)
allTransitionsPREV <- allTransitionsPREV[!(allTransitionsPREV$siteID=="UKFS"&allTransitionsPREV$year==2018),] #Not complete year
allTransitionsPREV <- allTransitionsPREV[!(allTransitionsPREV$siteID=="STEI"&allTransitionsPREV$year==2017),] #Not complete year

for(s in 1:length(site_names)){
  subDat <- allTransitionsPREV[allTransitionsPREV$siteID==site_names[s],]
  meanTrans$day15[s] <- mean(as.numeric(subDat$day15),na.rm=TRUE)
  meanTrans$sd15[s] <- sd(as.numeric(subDat$day15),na.rm=TRUE)
  meanTrans$day50[s] <- mean(as.numeric(subDat$day50),na.rm=TRUE)
  meanTrans$sd50[s] <- sd(as.numeric(subDat$day50),na.rm=TRUE)
  meanTrans$day85[s] <- mean(as.numeric(subDat$day85),na.rm=TRUE)
  meanTrans$sd85[s] <- sd(as.numeric(subDat$day85),na.rm=TRUE)
  
}
write.csv(meanTrans,file="historicalMeanTransitionDates.csv",row.names = FALSE)


