library("rjags")
library("ecoforecastR")
library("ncdf4")

source("randomWalkNullModelFunction.R")
###Random WalkNull Model Calculations
####Note: Currently this is not set up to run iteratively because I am not sure how the challenge is planning on doing this. 
####Note (continued): Hopefully someone who knows about how this will be done can use this code to do that

phenoDat <- read.csv("phenology-targets.csv.gz",header=TRUE)
sites <- unique(as.character(phenoDat$siteID))

forecast_length <- 35
predictions <- array(NA, dim = c(forecast_length, 1000, length(sites)))
parameters <- array(NA, dim = c(1000, length(sites)))

for(i in 1:length(sites)){
  
  message(paste0("forecasting site: ",sites[i]))
  
  forecast_length <- 35
  
  sitePhenoDat <- phenoDat[phenoDat$siteID==sites[i],]
  sitePhenoDat$time <- lubridate::as_date(sitePhenoDat$time)
  full_time <- tibble::tibble(time = seq(min(sitePhenoDat$time), max(sitePhenoDat$time) + lubridate::days(forecast_length), by = "1 day")) 
  forecast_start_index <- which(full_time$time == max(sitePhenoDat$time) + lubridate::days(1))
  d <- tibble::tibble(time = sitePhenoDat$time,
                      p=as.numeric(sitePhenoDat$gcc_90),
                      p.prec=1/((as.numeric(sitePhenoDat$gcc_sd))^2)) 
  d <- dplyr::full_join(d, full_time)
  #gap fill the missing precisions by assigning them the average sd for the site
  d$p.prec[!is.finite(d$p.prec)] <- NA
  d$p.prec[is.na(d$p.prec)] <- mean(d$p.prec,na.rm=TRUE)
  d$N <- length(d$p)
  data <- list(p = d$p,
               p.prec = d$p.prec,
               N = length( d$p))
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
  thinAmount <- round(nrow(out.mat)/1000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  
  predictions[, , i] <- as.matrix(out.burn$predict)[ ,forecast_start_index:length(d$time)]
  parameters[ ,i] <- as.matrix(out.burn$params)
  
}

team_name <- "EFInull"

forecast_time <- d$time[forecast_start_index:length(d$time)]

##Put in EFI standard form (based on EFI standards logistic-metadata-example vignette)
##Forecast Identifiers (please change to what is needed)

forecast_project_id <- team_name
forecast_model_id <- "v0.1"
forecast_iteration_id <- Sys.time()

#plot(predictions[1,], type = 'l', ylim = range(c(predictions)))
#for(i in 2:nrow(predictions)){
#  points(predictions[i,], type = "l")
#}

###Define Dimensions
timedim <- ncdim_def("time",
                     units = paste('days since', as.Date(forecast_time[1])),
                     vals = as.numeric(forecast_time - as.Date(forecast_time[1])),
                     longname='time')
ensdim <- ncdim_def("ensemble",
                    units="",
                    vals=1:nrow(as.matrix(out.burn$params)),
                    longname = "ensemble member")

sitedim <- ncdim_def("site",
                     units="",
                     vals=1:length(sites),
                     longname = "siteID")

dimnchar   <- ncdim_def("nchar",   "", 1:4, create_dimvar=FALSE )
## quick check that units are valid
udunits2::ud.is.parseable(timedim$units)
udunits2::ud.is.parseable(ensdim$units)
udunits2::ud.is.parseable(sitedim$units)

###Define Variables
def_list <- list()
def_list[[1]] <- ncvar_def(name = "gcc_90",
                           units = "",
                           dim = list(ensdim,timedim,sitedim),
                           longname = "90% quantile of daily green chromatic coordinate")
def_list[[2]] <- ncvar_def(name = "p.proc",
                           units = "",
                           dim = list(ensdim, sitedim),
                           longname = "Process precision parameter")
def_list[[3]] <- ncvar_def(name = "siteID",
                           units = "",
                           dim = list(dimnchar, sitedim),
                           longname = "siteID",
                           prec="char")

###Open netCDF file

forecast_file_name <- paste0("phenology-",lubridate::as_date(forecast_time[1] - lubridate::days(1)),"-",team_name,".nc")

ncout <- nc_create(forecast_file_name,def_list,force_v4=T)

###Fill in output data
ncvar_put(ncout,def_list[[1]], predictions) #Forecasted gcc_90
ncvar_put(ncout,def_list[[2]], parameters) #Forecasted parameter values
ncvar_put(ncout,def_list[[3]], sites) #Forecasted parameter values

## Global attributes (metadata)
ncatt_put(ncout,0,"forecast_project_id", as.character(forecast_project_id), 
          prec =  "text")
ncatt_put(ncout,0,"forecast_model_id",as.character(forecast_model_id), 
          prec =  "text")
ncatt_put(ncout,0,"forecast_iteration_id",as.character(forecast_iteration_id), 
          prec =  "text")
nc_close(ncout)   ## make sure to close the file
