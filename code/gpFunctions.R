# Internal functions of Gaussian Process main script

library(Matrix)
library(compiler)
#library(gputools)
# library(rpud)
# library(HiPLARM)

# Common functions

formatData <- function( u, y, lu, ly, flagPresent=FALSE ){
  n <- length(y)
  l <- max(lu, ly)
  if( l > n ){
    l <- n
  }
  
  x <- NULL
  if( flagPresent ){
    x <- matrix( 0, nrow=n-l, ncol=lu+ly+1 )  
  } else{
    x <- matrix( 0, nrow=n-l, ncol=lu+ly )  
  }
  
  out <- rep( 0, nrow(x) )
  for( i in (l+1):n ){
    if( flagPresent ){
      iniU <- i-1+1
    } else{
      iniU <- i-1  
    }    
    endU <- i-lu
    iniY <- i-1
    endY <- i-ly
    if( lu <= 0 ){
      x[i-l,] <- y[iniY:endY]  
    } else if( ly <= 0 ){
      x[i-l,] <- u[iniU:endU]    
    } else{
      x[i-l,] <- c( y[iniY:endY], u[iniU:endU] )    
    }      
    out[i-l] <- y[i]
  }  
  list( x=x, y=out ) 
}

normalizeData <- function( x, mean=NULL, std=NULL ){
  flagVector <- FALSE
  if( is.vector(x) || length(dim(x) == 1) ){
    x <- as.matrix(x)
    flagVector <- TRUE
  }
  if( is.null( mean ) ){
    mean <- colMeans(x)
  }
  if( is.null( std ) ){
    std <- sqrt( apply( x, 2, var ) )
  }  
  
  res <- t( apply( x, 1, function(d) (d - mean)/std ) )
  if( flagVector ){
    res <- as.vector(res)  
  }
  list( data=res, mean=mean, std=std )
}

gpOptimizeOrder2 <- function( u, y, lMax=5, kernel ){
  
  n <- length(u)
  trainRatio <- 0.7
  estimIndex <- 1:round(trainRatio*n)
  validateIndex <- (round(trainRatio*n)+1):n
  
  uEstim <- u[estimIndex]
  yEstim <- y[estimIndex]
  
  uValidate <- u[validateIndex]
  yValidate <- y[validateIndex]    
  
  grid <- expand.grid( 1:lMax, 1:lMax )
  grid <- grid[ which(grid[2] >= grid[1] ), ]
  bestLu <- 1; bestLy <- 1; bestValue <- Inf  
  
  for( i in 1:nrow(grid) ){    
    lu <- grid[i,1]; ly <- grid[i,2]
    print(paste("Avaliando", i, "/", nrow(grid), ", lu=", lu, ", ly=", ly))
    
    model <- gpRegression( u=uEstim, y=yEstim, lu=lu, ly=ly, kernel=kernel, optimizeOrder=FALSE, verbose=FALSE )
#     validation <- gpPredict( u=uValidate, y=yValidate, m=model, verbose=FALSE )               
#     err <- yValidate[(max(lu,ly)+1):length(uValidate)] - validation$est 
    validation <- gpSimulation( u=uValidate, y=yValidate, m=model, verbose=FALSE )               
    err <- yValidate - validation$est 
    
    mse <- mean(err^2)
    k <- length(model$par); n <- length(yValidate) 
    value <- 2*n*log(mse) + 2*k + (2*k*(k+1)) / (n-k-1)
    print(paste("AICc =", value))
    if( value < bestValue ){
      bestValue <- value
      bestLu <- lu; bestLy <- ly
    }    
  }  
  
  list( lu=bestLu, ly=bestLy )
}

gpOptimizeOrder <- function( u, y, lMax=5, kernel ){
    
  grid <- expand.grid( 1:lMax, 1:lMax )
  grid <- grid[ which(grid[2] >= grid[1] ), ]
  bestLu <- 1; bestLy <- 1; bestValue <- Inf  
  
  for( i in 1:nrow(grid) ){
    lu <- grid[i,1]; ly <- grid[i,2]
    print(paste("Avaliando", i, "/", nrow(grid), ", lu=", lu, ", ly=", ly))
    
    formatted <- formatData( u, y, lu, ly )
    xTrain <- formatted$x 
    yTrain <- formatted$y
    
    model <- gpRegression( x=xTrain, y=yTrain, lu=lu, ly=ly, kernel=kernel, optimizeOrder=FALSE, verbose=FALSE )
    
    L <- calcLogLikelihood( log(model$par), xTrain, yTrain, eps=eps )
    if( is.na(L) ){ next }    
    k <- length(model$par); n <- length(yTrain) 
    value <- -2*L + 2*k + (2*k*(k+1)) / (n-k-1)
    print(paste("AICc =", value))
    if( value < bestValue ){
      bestValue <- value
      bestLu <- lu; bestLy <- ly
    }    
  }  
  
  list( lu=bestLu, ly=bestLy )
}

gpSaveData <- function( data, name=NULL ){
  
  if( is.null(name) ){
    name <- format(Sys.time(), "%b %d %H:%M:%S %Y")
  }
  dir.create( file.path( getwd(), "gpSavedData" ), showWarnings=FALSE )
  save( data, file=paste( getwd(), "/gpSavedData/", name, ".dat", sep="" ) )
}

gpLoadData <- function( name, path=NULL ){
  
  if( is.null(path) ){
    load( file=paste( getwd(), "/gpSavedData/", name, ".dat", sep="" ) )  
  } else{
    load( file=paste( path, name, ".dat", sep="" ) )
  }
  
  if( exists("data") ){
    return( data)
  }
}

gpCheckData <- function( name, path=NULL ){
  
  if( is.null(path) ){
    filePath <- paste( getwd(), "/gpSavedData/", name, ".dat", sep="" )
  } else{
    filePath <- paste( path, name, ".dat", sep="" )
  }
  
  return( file.exists( filePath ) )
}

jitchol2 <- function( K, maxTries=5, text="" ){
  
  cholAux <- NA
  tryCatch( cholAux <- chol( K ), error=function(e) NA )
  if( any(is.na(cholAux)) == FALSE  ){ return(cholAux) }
  
  jit <- 10^-6 * mean(diag(K))
  identMat <- diag(NROW(K))
  for( i in 1:maxTries ){
    cholAux <- NA
    tryCatch( cholAux <- chol( K + jit*identMat ), error=function(e) NA )
    # if( any(is.na(cholAux)) == FALSE  ){ cat("JIT =", jit, text, "\n"); return(cholAux) }
    if( any(is.na(cholAux)) == FALSE  ){ return(cholAux) }
    jit <- 10*jit
  }
  # cat("FALHOU JIT =", jit, text, "\n")
  return(NA)
}

jitchol <- function( K, maxTries=5, text="" ){
  
  cholAux <- NA
  tryCatch( cholAux <- chol( K ), error=function(e) NA )
  if( any(is.na(cholAux)) == FALSE  ){ return(cholAux) }
  
  jit <- 10^-10
  identMat <- diag(nrow(K))
  for( i in 1:maxTries ){
    cholAux <- NA
    tryCatch( cholAux <- chol( K + jit*identMat ), error=function(e) NA )
    # if( any(is.na(cholAux)) == FALSE  ){ cat("JIT =", jit, text, "\n"); return(cholAux) }
    if( any(is.na(cholAux)) == FALSE  ){ return(cholAux) }
    jit <- 10*jit
  }
  # cat("JIT FAILED =", jit, text, "\n")
  return(NA)
}

