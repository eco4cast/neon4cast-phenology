##ESA 2021 Presentation Figures

##This script assumes that you have the following data files:
submittedForecasts <- read.csv("phenology_combined_forecasts.csv",header=TRUE)
allTransitions <- read.csv("allPhenologyTransitionData.csv",header=TRUE) #Created in calculatePhenoCamTransitionDates.R script

library(scales)
library(colorBlindness)
library(dplyr)

s=1
tms <- unique(submittedForecasts$team)
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")

#Figure 1: Dates of forecast and submission figure 

tranDates <- as.Date(allTransitions[,2],origin=as.Date("2020-12-31"))
tranDates <- c(tranDates,as.Date(allTransitions[,6],origin=as.Date("2020-12-31")))
challengeDays <- seq(as.Date("2021-02-01"),as.Date("2021-06-30"),"day")

jpeg(file="DatesOfForecastAndSubmissionFigure.jpg",width = 1000, height = 480, units = "px")#,height=3,width=10,units="inch",res=700)

par(mfrow=c(1,2),mai=c(1,2,0.3,0.1))

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

#Plot first team and then add on subsequent teams onto the graph
tm=1

plot(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20,
     xlab="Time",ylab="",ylim=c(0,length(tms)),bty="n",yaxt="n",main="Forecasted Days")
axis(side = 2,at=seq(1,length(tms)),labels=tms,pos=as.Date("2021-01-20"),las=1)

polygon(x=c(min(tranDates),min(tranDates),max(tranDates),max(tranDates)),y=c(-0.9,length(tms)+1,length(tms)+1,-0.9),col="chartreuse3",border=NA)
for(tm in 1:length(tms)){
  points(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20)
}

#Plot dates of submissions 
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

polygon(x=c(min(tranDates),min(tranDates),max(tranDates),max(tranDates)),y=c(-0.9,length(tms)+1,length(tms)+1,-0.9),col="chartreuse3",border=NA)
for(tm in 1:length(tms)){
  points(challengeDays[submittedForecastsYes[,tm]],rep(tm,length(challengeDays))[submittedForecastsYes[,tm]],pch=20)
}
#}
dev.off()

#Figure 2: Changes in forecasted values on transition dates
s=1

sitDat <- submittedForecasts[submittedForecasts$siteID==site_names[s],] #Subset by site
cls <- c("#004949","#000000",paletteMartin[3:15])

pdf(file="ForecastedValuesOnTransitionDates_presentationFigures.pdf",height=5,width=12)
par(mfrow=c(1,4),mai=c(0.5,0.5,1,0.1))

plotForecastedValuesOverTime <- function(highlightTms=NA){
  finalTms <- character()
  for(t in c(2,5,8)){ #Loops over the transition dates
    cl <- 1
    tranDate <- as.Date(allTransitions[s,t],origin=as.Date("2020-12-31"))
    vl <- allTransitions[s,(t+1)]
    sdVal <- allTransitions[s,(t+2)]
    vl <- rescale(vl,to=c(0,1),from=c(allTransitions$minimum[s],allTransitions$maximum[s])) ##Rescales gcc values between 0 and 1
    plot(x=numeric(),y=numeric(),type="l",xlim=c(-35,0),ylim=c(0,1),ylab="",xlab="Days Before Transition Date",main=paste(site_names[s],tranDate),bty="n")
    
    abline(h=vl,col="red",lwd=5,lty=2)
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
      if(nrow(organizedDat)>0){
        if(t==8){
          finalTms <- c(finalTms,as.character(tms[tm]))
        }
        sitTmMax <- max(tmSitDat$gcc_90[tmSitDat$statistic=="mean"],na.rm=TRUE)
        rescaledDat <- rescale(organizedDat$mean,to=c(0,1),from=c(allTransitions$minimum[s],sitTmMax)) #Rescales forecasted values between 0 and 1
        for(j in 1:length(rescaledDat)){ #Some values get scaled below 0 
          rescaledDat[j] <- max(rescaledDat[j],0)
        }
        if(is.na(highlightTms) || tm%in%highlightTms){ #No transparency if you do not want to highlight teams 
                                                       #or if the team is within the highlighted teams
          tF <- 1
          lwdVl <- 3
        }else{
          tF <- 0.2
          lwdVl <- 1
        }
        lines(as.Date(organizedDat$forecast_start_time)-tranDate,rescaledDat,col=scales::alpha(cls[cl],tF),lwd=lwdVl)
        cl <- cl + 1
      }
    }
  } 
  plot(x=numeric(),y=numeric(),type="l",xlim=c(-35,0),ylim=c(0,1),ylab="",xlab="",main="Legend",bty="n")
  legend("topleft",c("True Value",as.character(finalTms)),col=c("red",cls[1:length(finalTms)]),lty=c(2,rep(1,length(finalTms))),lwd=c(2,rep(3,length(finalTms))),bty = "n")
}
plotForecastedValuesOverTime()
plotForecastedValuesOverTime(highlightTms = 12)
plotForecastedValuesOverTime(highlightTms = 1)
plotForecastedValuesOverTime(highlightTms = 13)
plotForecastedValuesOverTime(highlightTms = c(2,11))
dev.off()
