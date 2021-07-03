
rm(list=ls())
set.seed(54)

## Load packages
library(RMediation)
library(dplyr)
library(boot)

## Define functions to generate data for simulation
generatedata <- function(apath, bpath, cPrimePath, nsamp){
  xdat = rnorm(nsamp)  # x variable
  mdat = apath*xdat + rnorm(nsamp)  # regression equation where x predicts m with error
  ydat = cPrimePath*xdat + bpath*mdat + rnorm(nsamp) # regression equation where x and m predict y with error
  data <- data.frame(xdat,mdat,ydat)
  return(data)
}

## Define function to run a linear regression (more efficient than R's built in lm() function)
fit.linear.regression <- function(outcome, predictors){    
  nsample <- NROW(predictors)
  Xmatrix <- cbind(rep(1, nsample),as.matrix(predictors))
  bmatrix <- solve(t(Xmatrix)%*%Xmatrix)%*%t(Xmatrix)%*%outcome
  yhat <- Xmatrix%*%bmatrix
  ematrix <- outcome-yhat
  ssresid <- t(ematrix)%*%ematrix
  msresid <- ssresid/(nsample-NROW(bmatrix))
  sigma <- as.numeric(msresid)*solve(t(Xmatrix)%*%Xmatrix)
  stderr <- sqrt(diag(sigma))
  stderr <- stderr[2]
  estimate <- bmatrix[2]
  tstat <- estimate/stderr
  p <- 2*(pt(abs(tstat),nsample-2,lower.tail=F))
  info <- (cbind(estimate,stderr,p))
  return(info)
}

regression.estimate <- function(outcome, predictors){    
  nsample <- NROW(predictors)
  Xmatrix <- cbind(rep(1, nsample),as.matrix(predictors))
  bmatrix <- solve(t(Xmatrix)%*%Xmatrix)%*%t(Xmatrix)%*%outcome
  estimate <- bmatrix[2]
  return(estimate)
}

## Define functions to run each of the tests of the indirect effect

causalSteps <- function(xdat, mdat, ydat, alpha){
  # calculate c path
  totalEffect <- fit.linear.regression(ydat, xdat)
  pc <- totalEffect[3]  #pc is the p-value for the total effect, step 1 in this method
  tao = totalEffect[1]
  # calculate a path
  XonM <- fit.linear.regression(mdat, xdat)
  pa <- XonM[3]
  # calculate b path
  MonY <- fit.linear.regression(ydat, cbind(mdat, xdat))
  pb <- MonY[3]
  XonYcontrollingM <- fit.linear.regression(ydat, cbind(xdat, mdat))
  taoPrime = XonYcontrollingM[1]
  if (pc > alpha){sigCS = 0
  } else if (pa > alpha){sigCS = 0
  } else if (pb > alpha){sigCS = 0
  } else if (taoPrime > tao){sigCS = 0
  } else {sigCS = 1}
  return(sigCS)
}

jointSignificance <- function(xdat, mdat, ydat, alpha){
  XonM <- fit.linear.regression(mdat, xdat)
  pa <- XonM[3]  #pa is the p-value for the effect of X on M, step 2 in this method
  if (pa < alpha){
    MonY <- fit.linear.regression(ydat, cbind(mdat, xdat))
    pb <- MonY[3] #pb is the p-value for the effect of M on Y controlling for X, step 3 in this method
    if (pb < alpha){
        sigJS = 1}
    else {sigJS = 0}}
  else {sigJS = 0}
  return(sigJS)
}

sobel <- function(xdat, mdat, ydat, alpha){
  # calculate a path
  XonM <- fit.linear.regression(mdat, xdat)
  estimatea <- XonM[1]
  stderra <- XonM[2]
  # calculate b path
  MonY <- fit.linear.regression(ydat, cbind(mdat, xdat))
  estimateb <- MonY[1]
  stderrb <- MonY[2]
  # calculate delta coefficient and test statistics
  delta <- sqrt(estimatea^2*stderrb^2+estimateb^2*stderra^2)
  tstatSobel <- (estimatea*estimateb)/delta
  pSobel <- 2*(1-pnorm(abs(tstatSobel)))
  # test significance
  if (pSobel < alpha){
    sigSobel = 1}
  else {sigSobel = 0}
  return(sigSobel)
}

prodclin <- function(xdat, mdat, ydat, alpha){
  # calculate a path
  XonM <- fit.linear.regression(mdat, xdat)
  estimatea <- XonM[1]
  stderra <- XonM[2]
  # calculate b path
  MonY <- fit.linear.regression(ydat, cbind(mdat, xdat))
  estimateb <- MonY[1]
  stderrb <- MonY[2]
  # use RMediation package to get a confidence interval
  interval <- medci(mu.x = estimatea, mu.y = estimateb, se.x = stderra, se.y = stderrb, alpha = alpha, type="prodclin")
  #If RMediation cannot calculate CI, output values that generated the error to .csv file and redo the iteration
  if (NA %in% interval$`97.5% CI`) {
    write.table(data.frame(estimatea, stderra, estimateb, stderrb), 
                "prodclinErrors.csv",
                append=TRUE, sep = ",",
                quote = FALSE,
                row.names=FALSE,
                col.names=FALSE)
    sigProdclin=999} 
  else{
  lowerbound <- interval$`97.5% CI`[1]
  upperbound <- interval$`97.5% CI`[2]
  if ((lowerbound > 0) | (upperbound < 0)){
    sigProdclin = 1}
  else {sigProdclin = 0}}
  return(sigProdclin)
}

percentileBootstrap <- function(data, alpha, nsamp, nbootstrap){
  # make a vector for the indirect effect estimates
  indirectEffectBootstrap <- rep(999,nbootstrap)
  for (e in 1:nbootstrap){
    # bootstrap the data
    bootstrapScores = sample_n(data, size = nsamp, replace = T)
    xBootstrap = bootstrapScores[,1]
    mBootstrap = bootstrapScores[,2]
    yBootstrap = bootstrapScores[,3]
    # calculate the bootstrapped a path
    estimateaBootstrap <- regression.estimate(mBootstrap,xBootstrap)
    # calculate the bootstrapped b path
    estimatebBootstrap <-regression.estimate(yBootstrap, cbind(mBootstrap, xBootstrap))
    # calculate the indirect effect
    indirectEffectBootstrap[e] = estimateaBootstrap*estimatebBootstrap}
  # create a confidence interval for the bootstrapped indirect effect
  indirectEffectBootstrap <- sort(indirectEffectBootstrap)
  lowerPB <- indirectEffectBootstrap[floor((alpha/2)*nbootstrap)]
  upperPB <- indirectEffectBootstrap[ceiling((1-alpha/2)*nbootstrap)]
  if ((lowerPB > 0) | (upperPB < 0)){
    sigPB = 1}
  else {sigPB = 0}
  return(sigPB)
}

BiasCorrectedBootstrap <- function(data, alpha, nsamp, nbootstrap){
  # make a vector for the indirect effect estimates
  indirectEffectBootstrap <- rep(999,nbootstrap)
  for (d in 1:nbootstrap){
    # bootstrap the data
    bootstrapScores = sample_n(data, size = nsamp, replace = T)
    xBootstrap = bootstrapScores[,1]
    mBootstrap = bootstrapScores[,2]
    yBootstrap = bootstrapScores[,3]
    # calculate the bootstrapped a path
    estimateaBootstrap <- regression.estimate(mBootstrap,xBootstrap)
    # calculate the bootstrapped b path
    estimatebBootstrap <-regression.estimate(yBootstrap, cbind(mBootstrap, xBootstrap))
    # calculate the indirect effect
    indirectEffectBootstrap[d] = estimateaBootstrap*estimatebBootstrap}
  ## calculate original indirect effect
  # calculate a path
  estimatea <- regression.estimate(mdat, xdat)
  # calculate b path
  estimateb <- regression.estimate(ydat, cbind(mdat,xdat))
  # original indirect effect
  indirectEffect = estimatea*estimateb
  # create a confidence interval for the bootstrapped indirect effect
  indirectEffectBootstrap <- sort(indirectEffectBootstrap)
  # bias-correction
  propUnderOriginal <- sum(indirectEffectBootstrap<indirectEffect)/nbootstrap
  z0 <- qnorm(propUnderOriginal)
  zcritl <- qnorm(alpha/2)
  zcritu <- qnorm(1-alpha/2)
  zl <- pnorm(2*z0+zcritl)
  zu <- pnorm(2*z0+zcritu)

  lowerBCB <- indirectEffectBootstrap[floor((zl)*nbootstrap)]
  upperBCB <- indirectEffectBootstrap[ceiling((zu)*nbootstrap)]
  if ((lowerBCB > 0) | (upperBCB < 0)) {sigBCB = 1}
  else {sigBCB=0}
  return(sigBCB)
}








## Set Parameters
apath = c(.14, .26, .39, .59)
bpath = c(.14, .26, .39, .59)
cPrimePath = c(.14, .39, .59, 0)
alpha = .05
desiredPower = .8
nbootstrap = 2000

samplesizesMatrix <- data.frame(matrix(ncol = 5, nrow = 4*4*4*6))  # 4 effect sizes on each a, b, and c, then 6 tests
row = 1


for (b in 1:length(bpath)){  # change to 1:length(bpath)
  for (a in 1:length(apath)){   
    
    
    #### Simulation
    for (cprime in 3:length(cPrimePath)){
      
      print(apath[a])
      print(bpath[b])
      print(cPrimePath[cprime])
      
      for (test in 1:6){
        
        print(test)
        
        whileCounter = 0
        powerEstimate = .5
        
        ## give a starting sample size estimate
        if (test == 1 & cPrimePath[cprime] == 0 & apath[a]==.14 & bpath[b]==.14){
          nsamp = 20000} 
        else if (apath[a]==.14 | bpath[b]==.14){
          nsamp = 300}
        else {
          nsamp = 50}
        
        ## set number of sims according to article
        if (test <= 4){
          nsims = 100000
          errorMargin = .001}
        else {
          nsims = 1000
          errorMargin = .005}
        
        if (apath[a]==.59 | bpath[b]==.59){  ## this is because there is no sample size that gets exactly .8 power at the largest effect size
          errorMargin = .01
        }
        
        while (powerEstimate < desiredPower-errorMargin | powerEstimate > desiredPower+errorMargin){
       
          # initialize counter variables
          sumCSrejection = 0
          sumJSrejection = 0
          sumSobelRejection = 0
          sumProdclinRejection = 0
          sumPBrejection = 0
          sumBCBrejection = 0

          for (i in 1:nsims){
            ## generate data
            data <- generatedata(apath[a], bpath[b], cPrimePath[cprime], nsamp)
            xdat <- data[,1]
            mdat <- data[,2]
            ydat <- data[,3]
    
            ## run tests
            # Causal Steps
            if (test == 1){
              CSsig <- causalSteps(xdat, mdat, ydat, alpha)
              sumCSrejection = sumCSrejection + CSsig}
            # Joint Significance
            if (test == 2){
              JSsig <- jointSignificance(xdat, mdat, ydat, alpha)
              sumJSrejection = sumJSrejection + JSsig}
            # Sobel
            if (test == 3){
              sobelSig <- sobel(xdat, mdat, ydat, alpha)
              sumSobelRejection = sumSobelRejection + sobelSig}
            # Prodclin
            if (test == 4){
              prodclinSig=999
              while(prodclinSig==999){
                data <- generatedata(apath[a], bpath[b], cPrimePath[cprime], nsamp)
                xdat <- data[,1]
                mdat <- data[,2]
                ydat <- data[,3]
              prodclinSig <- prodclin(xdat, mdat, ydat, alpha)}
              sumProdclinRejection = sumProdclinRejection + prodclinSig}
            # Percentile Bootstrap
            if (test == 5){
              pbSig <- percentileBootstrap(data, alpha, nsamp, nbootstrap)
              sumPBrejection = sumPBrejection + pbSig}
            # Bias-corrected Bootstrap
            if (test == 6){
              bcbSig <- BiasCorrectedBootstrap(data, alpha, nsamp, nbootstrap)
              sumBCBrejection = sumBCBrejection + bcbSig}
          }
          
          # save previous estimates for slope calculation later
          if (whileCounter > 0){
            previousPowerEstimate = powerEstimate}
    
          # calculate power for test
          if (test == 1){
            powerCS = sumCSrejection/nsims
            powerEstimate = powerCS}
          if (test == 2){
            powerJS = sumJSrejection/nsims
            powerEstimate = powerJS}
          if (test == 3){
            powerSobel = sumSobelRejection/nsims
            powerEstimate = powerSobel}
          if (test == 4){
            powerProdclin = sumProdclinRejection/nsims
            powerEstimate = powerProdclin}
          if (test == 5){
            powerPB = sumPBrejection/nsims
            powerEstimate = powerPB}
          if (test == 6){
            powerBCB = sumBCBrejection/nsims
            powerEstimate = powerBCB} 
          
          # sim will crash if new power estimate is the same as previous one
          if (whileCounter > 0){
            if (powerEstimate == previousPowerEstimate){powerEstimate = powerEstimate-.01}}
    
          
          if (whileCounter == 0){
            slope = (powerEstimate-0)/(nsamp-0)}
          else{slope = (powerEstimate-previousPowerEstimate)/(nsamp-previousNsamp)}
          
          previousNsamp = nsamp
          
          nsamp = abs(nsamp + ceiling((desiredPower-powerEstimate)/slope))
          
          # sim will crash if new sample size is the same as previous one
          if (nsamp == previousNsamp){nsamp = nsamp+1}
          
          print(previousNsamp)
          print(nsamp)
          print(powerEstimate)
        
          
          whileCounter = whileCounter+1
        }
          
        # save the sample size that achieves the desired level of power for the appropriate test
        samplesizesMatrix[row,] <- c(previousNsamp, test, apath[a], bpath[b], cPrimePath[cprime])
        write.csv(samplesizesMatrix, "samplesizes.csv")
        row = row + 1
      
      }
     }
    
  }
}


## Put together table
library(tidyverse)
dat <- read.csv("samplesizes.csv")
colnames(dat) <- c("row","n","test","apath","bpath","cpath")
dat <- dat %>% filter(!is.na(n))

table3 <- rbind(
  ## Causal Steps tao=0 row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$cpath==0 & dat$test==1] %>% mean() %>% round(0)),
  ## Causal Steps tao=.14 row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$cpath==.14 & dat$test==1] %>% mean() %>% round(0)),
  ## Causal Steps tao=.39 row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$cpath==.39 & dat$test==1] %>% mean() %>% round(0)),
  ## Causal Steps tao=.59 row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$cpath==.59 & dat$test==1] %>% mean() %>% round(0)),
  ## joint significance test row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$test==2] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$test==2] %>% mean() %>% round(0)),
  ## Sobel test row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$test==3] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$test==3] %>% mean() %>% round(0)),
  ## PRODCLIN row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$test==4] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$test==4] %>% mean() %>% round(0)),
  ## Percentile Bootstrap row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$test==5] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$test==5] %>% mean() %>% round(0)),
  ## BCB row
  c(dat$n[dat$apath==.14 & dat$bpath==.14 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.26 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.39 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.14 & dat$bpath==.59 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.14 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.26 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.39 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.26 & dat$bpath==.59 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.14 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.26 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.39 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.39 & dat$bpath==.59 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.14 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.26 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.39 & dat$test==6] %>% mean() %>% round(0),
    dat$n[dat$apath==.59 & dat$bpath==.59 & dat$test==6] %>% mean() %>% round(0)))
View(table3)
write.table(table3, 
            "OutputTable.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=FALSE)



