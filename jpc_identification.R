closeAllConnections()
rm(list=ls())
library(compiler)
library(ggplot2)
source("data/loadData.R")
source("code/validation.R")
source("code/figuresFunctions.R")

###############################################################################

# Change the current path to the folder containing this file and run the scriot.

### Choose the options below by commenting and uncommenting the related lines.

### Choose model
# algorithName <- "gp"; source("code/gpnarx.R") # GP-NARX
algorithName <- "gp"; source("code/rgpt.R") # RGP-t

### Choose dataset
# dataName <- "identificationExample" # Artificial example without outliers
dataName <- "identificationExample_30outt" # Artificial example with 30% of outliers sampled from a Student-t distribution


###############################################################################

### Do not change the settings below

trainRatio <- 0.5
normalizeSeries <- TRUE
# Number of hidden layers for RGP-t
numHiddenLayers <- 1
# Number of pseudo-inputs for RGP-t
M <- 15

m <- NULL
identifySystem <- function( algorithName, dataName ){
  
  set.seed(123)
  
  # Load data
  data <- loadData( dataName )
  dataName <- data$dataName
  u <- data$u
  y <- data$y
  uNew <- data$uNew
  yNew <- data$yNew
  lu <- data$lu
  ly <- data$ly
  
  if( is.null(data$continueSimulation) == FALSE ){
    continueSimulation <<- data$continueSimulation
  }
  
  if( is.null(yNew) ){
    continueSimulation <<- TRUE
    uTotal <- u
    yTotal <- y
    n <- length(uTotal)
    estimIndex <- 1:round(trainRatio*n)
    avalIndex <- (round(trainRatio*n)+1):n
    
    u <- uTotal[estimIndex]
    y <- yTotal[estimIndex]
    
    uNew <- uTotal[avalIndex]
    yNew <- yTotal[avalIndex]
  }
  
  if( normalizeSeries ){
    uMean <- mean(u); uSd <- sd(u)
    yMean <- mean(y); ySd <- sd(y)
    uOriginal <- u; yOriginal <- y
    u <- u - uMean; u <- u/uSd
    y <- y - yMean; y <- y/ySd
    uNewOriginal <- uNew; yNewOriginal <- yNew
    uNew <- uNew - uMean; uNew <- uNew/uSd
    yNew <- yNew - yMean; yNew <- yNew/ySd
  } else{
    uOriginal <- u; yOriginal <- y
    uNewOriginal <- uNew; yNewOriginal <- yNew  
  }
  
  if( is.null(lu) ){
    lu <- 2
    ly <- 2
  }
  l <- max(lu,ly)
  
  algorithName <- toupper(algorithName)
  regressionFunction <- gpRegression
  predictionFunction <- gpPredict
  simulationFunction <- gpSimulation
  
  assign( "uTrain", u, envir=.GlobalEnv )
  assign( "yTrain", y, envir=.GlobalEnv )
  assign( "uNew", uNew, envir=.GlobalEnv )
  assign( "yNew", yNew, envir=.GlobalEnv )
  
  print("", quote=FALSE)
  print(paste("Model estimation", algorithName, "..."), quote=FALSE)  
  m <- regressionFunction( u=u, y=y, lu=lu, ly=ly, par=NULL, M=M, numHiddenLayers=numHiddenLayers, robustFlag=TRUE, verbose=TRUE ) 
  m <<- m
  
  print("", quote=FALSE)
  print("Free simulation with test data...", quote=FALSE)
  simul <- NA
  tryCatch( simul <- simulationFunction( m=m, u=uNew, y=yNew, continueSimulation=continueSimulation ), error=function(e) NA )
  if( any(is.na(simul)) == FALSE ){
    if( normalizeSeries ){ 
      simul$est <- simul$est*ySd + yMean
      if( ! is.null( simul$var ) ){ simul$var <- simul$var*(ySd^2) } 
    }
    estSimul <- simul$est#[(l+1):length(simul$est)]
    if( is.null(simul$var) == FALSE ){
      varSimul <- simul$var#[(l+1):length(simul$var)]
      varSimul <- ifelse( is.na(varSimul) | varSimul < 0, min(varSimul[varSimul>0],na.rm=TRUE), varSimul )
    } else{
      varSimul <- NULL  
    }
    
    if( length(estSimul) < length(yNewOriginal) ){
      yNewOriginalComp <- yNewOriginal[(l+1):length(simul$est)]
    } else{
      yNewOriginalComp <- yNewOriginal
    }
    metricsSimulation <- calculateMetrics( yNewOriginalComp, estSimul, varSimul )  
    
    plotValidation( yNewOriginalComp, estSimul, varSimul, label="Simulation",
                    figName=paste(algorithName, "-", dataName,"- Simulation"),
                    sub=paste("RMSE =", format(round(metricsSimulation$rmse, 4), nsmall = 4),
                              ifelse(algorithName=="GP", paste(", NLD =", format(round(metricsSimulation$ld, 4), nsmall = 4)),"")) )
    
    res <- list( m=m, estSimul=estSimul, varSimul=varSimul, yNew=yNewOriginal, metricsSimulation=metricsSimulation )
  }
  
  return( res )
}

res <- identifySystem( algorithName, dataName )
