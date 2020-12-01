###Random Walk Null Model
##' Creates a random walk phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, and data$n
##' @param nchain The desired number of chains in the MCMC
##' @param priorCal If a calibration period has been performed enter the priors here in the form of a list (if not do not include)
##' @export
##' @import rjags
##' @import coda
randomWalkPhenoModel <- function(data,nchain,priorCal=FALSE){
  ##Set priors
  if(typeof(priorCal)==typeof(FALSE)){ ##Done when there was not a calibration period performed separately (or if this is the calibration)
    data$s1.proc <- 1262.626
    data$s2.proc <- 50.50505
    data$x1.a <- 1 #Done to keep distribution close to 0 (over 75% of the data <0.05)
    data$x1.b <- 30
    #data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    #data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
  }
  
  ###JAGS model
  RandomWalk = "
  model{
  
  #### Data Models
  for(i in 1:N){
    p[i] ~ dnorm(x[i],p.prec[i])
  }
  
  #### Process Model
  for(i in 2:N){
    xl[i]~dnorm(x[i-1],p.proc)
    x[i] <- max(0, min(1,xl[i]))
  }
  
  #### Priors
  x[1] ~ dbeta(x1.a,x1.b)
  #p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  }
  "
  
  ###Create the JAGS model using the basic RandomWalk Model
  
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           n.chains = nchain)
  return(j.model)
}

