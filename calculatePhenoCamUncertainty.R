##' Calculate the uncertainty (standard deviation) on daily PhenoCam gcc_90 by bootstrap
##'
##' @param dat PhenoCam data dataframe from roistats file
##' @param dates Vector of desired dates to calculate standard deviation for
##' @export
calculate.phenocam.uncertainty <- function(dat,dates) {
  sds <- array(NA,dim = c(length(dates), 2))
  nboot <- 50
  for(d in 1:length(dates)){
    dailyDat <- dat[dat$date==dates[d],]
    if(nrow(dailyDat)>0){
      dailyDat <- dailyDat[!is.na(dailyDat$gcc),]
      nrows <- nrow(dailyDat)
      gcc_90s <- rep(NA,nboot)
      rcc_90s <- rep(NA,nboot)
      for(j in 1:nboot){
        gcc_90s[j] <- quantile(dailyDat$gcc[sample(x = 1:nrows,size = nrows,replace = T)],0.90)
        rcc_90s[j] <- quantile(dailyDat$rcc[sample(x = 1:nrows,size = nrows,replace = T)],0.90)
      }
      sds[d,1] <- sd(gcc_90s)
      sds[d,2] <- sd(rcc_90s)
    }else{
      sds[d, 1] <- NA
      sds[d, 2] <- NA
    }
  }
  return(sds)
}