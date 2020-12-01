library("rjags")
library("ecoforecastR")
library("ncdf4")

source("randomWalkNullModelFunction.R")
###Random WalkNull Model Calculations
####Note: Currently this is not set up to run iteratively because I am not sure how the challenge is planning on doing this. 
####Note (continued): Hopefully someone who knows about how this will be done can use this code to do that

phenoDat <- read.csv("phenology-targets.csv.gz",header=TRUE)
sites <- unique(as.character(phenoDat$siteID))
for(i in 1:length(sites)){
  sitePhenoDat <- phenoDat[phenoDat$siteID==sites[i],]
  data <- list(p=as.numeric(sitePhenoDat$gcc_90),p.prec=1/((as.numeric(sitePhenoDat$gcc_sd))**2)) 
  #gap fill the missing precisions by assigning them the average sd for the site
  data$p.prec[!is.finite(data$p.prec)] <- NA
  data$p.prec[is.na(data$p.prec)] <- mean(data$p.prec,na.rm=TRUE)
  data$N <- length(data$p)
  j.model <- randomWalkPhenoModel(data=data,nchain=5,priorCal=FALSE)
  variableNames <- c("x","p.proc") #x is the latent variable of gcc_90 and p.proc is the process precision 
  jags.out   <- coda.samples (model = j.model,
                              variable.names = variableNames,
                              n.iter = 20000)
  
  #Split output into parameters and state variables and calculat/remove burnin 
  out = list(params=NULL,predict=NULL) 
  mfit = as.matrix(jags.out,chains=TRUE)
  pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
  chain.col = which(colnames(mfit)=="CHAIN")
  out$params = ecoforecastR::mat2mcmc.list(mfit[,-pred.cols])
  GBR <- gelman.plot(out$params)
  burnin <- GBR$last.iter[tail(which(any(GBR$shrink[,,2]>1.05,1)),1)+1]

  var.burn <- window(jags.out,start=burnin)
  out.burn = list(params=NULL,predict=NULL)
  mfit = as.matrix(var.burn,chains=TRUE)
  pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
  chain.col = which(colnames(mfit)=="CHAIN")
  out.burn$params = ecoforecastR::mat2mcmc.list(mfit[,-pred.cols])
  out.burn$predict = ecoforecastR::mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
  
  ##Thin (not sure what the standards think about this)
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/5000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  
  ##Put in EFI standard form (based on EFI standards logistic-metadata-example vignette)
  ##Forecast Identifiers (please change to what is needed)
  
  forecast_project_id <- paste("PhenoForecastRandomWalk_",sitePhenoDat$siteID[1],sep="")
  forecast_model_id <- "v0.1"
  forecast_iteration_id <- Sys.time()
  times <- as.Date(sitePhenoDat$time)
  
  ###Define Dimensions
  timedim <- ncdim_def("time",
                       units = paste('days since', as.Date(times[1])),
                       vals = as.numeric(times - as.Date(times[1])),
                       longname='time')
  ensdim <- ncdim_def("ensemble",
                      units="",
                      vals=1:nrow(as.matrix(out.burn$params)),
                      longname = "ensemble member")
  ## quick check that units are valid
  udunits2::ud.is.parseable(timedim$units)
  udunits2::ud.is.parseable(ensdim$units)
  
  ###Define Variables
  def_list <- list()
  def_list[[1]] <- ncvar_def(name = "gcc_90",
                             units = "",
                             dim = list(ensdim,timedim),
                             longname = "90% quantile of daily green chromatic coordinate")
  def_list[[2]] <- ncvar_def(name = "p.proc",
                             units = "",
                             dim = list(ensdim),
                             longname = "Process precision parameter")

  ###Open netCDF file
  ncout <- nc_create(paste("pheno-forecast-random-walk_",sitePhenoDat$siteID[1],".nc",sep=""),def_list,force_v4=T)
  
  ###Fill in output data
  ncvar_put(ncout,def_list[[1]], as.matrix(out.burn$predict)) #Forecasted gcc_90
  ncvar_put(ncout,def_list[[2]], as.matrix(out.burn$params)) #Forecasted parameter values
  
  ## Global attributes (metadata)
  ncatt_put(ncout,0,"forecast_project_id", as.character(forecast_project_id), 
            prec =  "text")
  ncatt_put(ncout,0,"forecast_model_id",as.character(forecast_model_id), 
            prec =  "text")
  ncatt_put(ncout,0,"forecast_iteration_id",as.character(forecast_iteration_id), 
            prec =  "text")
  nc_close(ncout)   ## make sure to close the file
}