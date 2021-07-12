##Figure Ideas:

#Something about showing when people submitted forecasts.
submittedForecasts <- read.csv("phenology_combined_forecasts.csv",header=TRUE)

tms <- unique(submittedForecasts$team)
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
allTransitions <- read.csv("allPhenologyTransitionData.csv",header=TRUE)
tranDates <- as.Date(allTransitions[,2],origin=as.Date("2020-12-31"))
tranDates <- c(tranDates,as.Date(allTransitions[,6],origin=as.Date("2020-12-31")))
challengeDays <- seq(as.Date("2021-02-01"),as.Date("2021-06-30"),"day")

jpeg(file="DatesOfForecastAndSubmissionFigure.jpg",width = 1000, height = 480, units = "px")#,height=3,width=10,units="inch",res=700)

par(mfrow=c(1,2))
par(mai=c(1,2,0.3,0.1))
s <- 1
#for(s in 1:length(site_names)){
  submittedForecastsYes <- data.frame(matrix(ncol=length(tms),nrow=length(challengeDays)))
  for(tm in 1:length(tms)){
    tmDat <- submittedForecasts[submittedForecasts$team==tms[tm],] #Subset by team
    tmSitDat <- tmDat[tmDat$siteID==site_names[s],] #Subset by site
    uniqueTimes <- unique(as.Date(tmSitDat$time))
    for(d in 1:length(challengeDays)){
      submittedForecastsYes[d,tm] <- challengeDays[d] %in% uniqueTimes
    }
  }
  colnames(submittedForecastsYes) <- tms
  
  tm=1
  
  plot(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20,
       xlab="Time",ylab="",ylim=c(0,length(tms)),bty="n",yaxt="n",main="Forecasted Days")
  axis(side = 2,at=seq(1,length(tms)),labels=tms,pos=as.Date("2021-01-20"),las=1)
  
  polygon(x=c(min(tranDates),min(tranDates),max(tranDates),max(tranDates)),y=c(-0.9,length(tms)+1,length(tms)+1,-0.9),col="chartreuse3",border=NA)
  for(tm in 1:length(tms)){
    points(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20)
  }
  
  submittedForecastsYes <- data.frame(matrix(ncol=length(tms),nrow=length(challengeDays)))
  
  for(tm in 1:length(tms)){
    tmDat <- submittedForecasts[submittedForecasts$team==tms[tm],] #Subset by team
    tmSitDat <- tmDat[tmDat$siteID==site_names[s],] #Subset by site
    uniqueTimes <- unique(as.Date(tmSitDat$forecast_start_time))
    for(d in 1:length(challengeDays)){
      submittedForecastsYes[d,tm] <- challengeDays[d] %in% uniqueTimes
    }
  }
  colnames(submittedForecastsYes) <- tms
  
  tm=1
  par(mai=c(1,0.1,0.3,2))
  plot(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20,
       xlab="Time",ylab="",ylim=c(0,length(tms)),bty="n",yaxt="n",main="Submission Days",xlim=range(challengeDays))
  #axis(side = 2,at=seq(1,length(tms)),labels=tms,pos=as.Date("2021-01-20"),las=1)
  
  polygon(x=c(min(tranDates),min(tranDates),max(tranDates),max(tranDates)),y=c(-0.9,length(tms)+1,length(tms)+1,-0.9),col="chartreuse3",border=NA)
  for(tm in 1:length(tms)){
    points(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20)
  }
#}
dev.off()

##Need a figure about transition dates

library(tidyverse)
library(tools)

c24 <- c(
  "dodgerblue2",
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

submittedForecasts <- read.csv("phenology_combined_forecasts.csv",header=TRUE)
allTransitions <- read.csv("allPhenologyTransitionData.csv",header=TRUE) #Created in calculatePhenoCamTransitionDates.R script

tms <- unique(submittedForecasts$team)
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")

pdf(file="ForecastedValuesOnTransitionDates.pdf",height=10,width=15)
par(mfrow=c(1,3))

for(s in 1:length(site_names)){
  print(site_names[s])
  sitDat <- submittedForecasts[submittedForecasts$siteID==site_names[s],] #Subset by site
  for(t in c(2,5,8)){ #Loops over the transition dates
    tranDate <- as.Date(allTransitions[s,t],origin=as.Date("2020-12-31"))
    vl <- allTransitions[s,(t+1)]
    sdVal <- allTransitions[s,(t+2)]
    plot(x=tranDate,NA,type="l",xlim=c(tranDate-35,tranDate),ylim=c(0.3,0.5),ylab="GCC",xlab="Time",main=paste(site_names[s],tranDate))
    if(t==2){
      legend("topleft",as.character(tms),col=c24[1:length(tms)],lty=rep(1,length(tms)),lwd=rep(2,length(tms)))
    }
    polygon(c(tranDate-40,tranDate-40,tranDate+5,tranDate+5),c((vl-1.96*sdVal),(vl+1.96*sdVal),(vl+1.96*sdVal),(vl-1.96*sdVal)),col=alpha("red",0.5),border=FALSE)
    abline(h=vl,col="red",lwd=2)
    for(tm in 1:length(tms)){
      tmSitDat <- sitDat[sitDat$team==tms[tm],] #Subset by team
      
      predDat <- tmSitDat[as.Date(tmSitDat$time)==tranDate,] #Subset of the forecasts that forecasted the transition date 
      
      #Organizing data to put upper95, lower95, and mean on the same row 
      organizedDat <- predDat[predDat$statistic=="mean",c("time","forecast_start_time","gcc_90")]
      organizedDat <- organizedDat %>% rename("mean"="gcc_90")
      lowDat <- predDat[predDat$statistic=="lower95",c("time","forecast_start_time","gcc_90")]
      lowDat <- lowDat %>% rename("lower95"="gcc_90")
      highDat <- predDat[predDat$statistic=="upper95",c("time","forecast_start_time","gcc_90")]
      highDat <- highDat %>% rename("upper95"="gcc_90")
      organizedDat <- merge(organizedDat,lowDat,by=c("time","forecast_start_time"))
      organizedDat <- merge(organizedDat,highDat,by=c("time","forecast_start_time"))
      lines(as.Date(organizedDat$forecast_start_time),organizedDat$mean,col=c24[tm])
    }
  }
}

dev.off()

