##Script to evaluate how far in advance different forecasts were able to estimate the GCC of transition dates within various
##levels of accuracy. 

submittedForecasts <- read.csv("phenology_combined_forecasts.csv",header=TRUE)
allTransitions <- read.csv("allPhenologyTransitionData.csv",header=TRUE) #Created in calculatePhenoCamTransitionDates.R script

tms <- unique(submittedForecasts$team)
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
cutOffs <- c(0.10,0.15,0.20,0.3)
limitNames <- c("r10","r15","r20","r30","UQ") #Names to indicate that the team's predicted value fell within 10%, 15%, 20%, 30% of the fitted model; UQ indicates it fell within the submitted model's 95% CI
transNames <- c("SOS","MOS","EOS") #Start of season; middle of season; end of season 

allSuccessfulDates <- data.frame(matrix(ncol=5,nrow=0))
colnames(allSuccessfulDates) <- c("team","siteID","transition","limit","dayInAdvance")
allSuccessfulDates$team <- as.character(allSuccessfulDates$team)
allSuccessfulDates$siteID <- as.character(allSuccessfulDates$siteID)
allSuccessfulDates$dayInAdvance <- as.numeric(allSuccessfulDates$dayInAdvance)
allSuccessfulDates$transition <- as.character(allSuccessfulDates$transition)
allSuccessfulDates$limit <- as.character(allSuccessfulDates$limit)

for(tm in 1:length(tms)){
  tmDat <- submittedForecasts[submittedForecasts$team==tms[tm],] #Subset by team
  for(s in 1:length(site_names)){
    tmSitDat <- tmDat[tmDat$siteID==site_names[s],] #Subset by site
    for(t in c(2,4,6)){ #Loops over the transition dates
      tranDate <- as.Date(allTransitions[s,t],origin=as.Date("2020-12-31"))
      vl <- allTransitions[s,(t+1)]
      
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
      
      if(nrow(organizedDat)>0){ #check to see if that team made a forecast for this date (and if not, not add any row(s) for that team for the transition date)
        for(c in 1:length(cutOffs)){
          successRows <- organizedDat[((cutOffs[c]*allTransitions$rangeLength[s]+vl) >organizedDat$mean & 
                                         (vl - cutOffs[c]*allTransitions$rangeLength[s])<organizedDat$mean),]
          if(nrow(successRows)>1){
            teamSuccessfulDates <- tranDate - as.Date(successRows$forecast_start_time) #Calculate how many days before the transition date, the forecast was made
            numDates <- length(teamSuccessfulDates)
            newDat <- cbind(rep(as.character(tms[tm]),numDates),rep(as.character(site_names[s]),numDates),
                            rep(transNames[t/2],numDates),
                            rep(as.character(limitNames[c]),numDates),as.numeric(teamSuccessfulDates))
            
            colnames(newDat) <- c("team","siteID","transition","limit","dayInAdvance")
            allSuccessfulDates <- rbind(allSuccessfulDates,newDat)
            allSuccessfulDates$dayInAdvance <- as.numeric(allSuccessfulDates$dayInAdvance)
            
          }else{
            allSuccessfulDates <- rbind(allSuccessfulDates,c(as.character(tms[tm]),site_names[s],transNames[t/2],
                                                             limitNames[c],0))
            allSuccessfulDates$team <- as.character(allSuccessfulDates$team) ##R kept changing classes to factors so these are to maintain the classes as not factors
            allSuccessfulDates$dayInAdvance <- as.numeric(allSuccessfulDates$dayInAdvance)
            allSuccessfulDates$transition <- as.character(allSuccessfulDates$transition)
            allSuccessfulDates$limit <- as.character(allSuccessfulDates$limit)
            allSuccessfulDates$siteID <- as.character(allSuccessfulDates$siteID)
          }
        }
        successRows <- organizedDat[(vl>organizedDat$lower95 && vl<organizedDat$upper95),]
        if(nrow(successRows)>0){
          teamSuccessfulDates <- tranDate - as.Date(successRows$forecast_start_time)
          numDates <- length(teamSuccessfulDates)
          newDat <- cbind(rep(as.character(tms[tm]),numDates),rep(as.character(site_names[s]),numDates),
                          rep(transNames[t/2],numDates),
                          rep(as.character(limitNames[5]),numDates),as.numeric(teamSuccessfulDates))
          
          colnames(newDat) <- c("team","siteID","transition","limit","dayInAdvance")
          allSuccessfulDates <- rbind(allSuccessfulDates,newDat)
          
        }else{
          allSuccessfulDates$team <- as.character(allSuccessfulDates$team)
          allSuccessfulDates$dayInAdvance <- as.numeric(allSuccessfulDates$dayInAdvance)
          allSuccessfulDates$transition <- as.character(allSuccessfulDates$transition)
          allSuccessfulDates$limit <- as.character(allSuccessfulDates$limit)
          allSuccessfulDates$siteID <- as.character(allSuccessfulDates$siteID)
          allSuccessfulDates <- rbind(allSuccessfulDates,c(as.character(tms[tm]),site_names[s],transNames[t/2],
                                                           limitNames[5],0))
        }
      }
    }
  }
}

write.csv(allSuccessfulDates,file="teamSuccessfulTransitionDates.csv",row.names=FALSE)
