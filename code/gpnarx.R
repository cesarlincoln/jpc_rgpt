library(Matrix)
library(compiler)
source('code/gpFunctions.R', chdir=TRUE)
source('code/gpKernels.R', chdir=TRUE)
source('code/gpLik.R', chdir=TRUE)

eps <-  10^-6 # 10^-10
generateRandomParams <- NULL
calcKernelMatrix <- NULL
calcGradLogLikelihood <- NULL
calcRegularizationPenalty <- NULL

calcKernelMatrix <- cmpfun( function( x1, x2, par=NULL ){  
  if( is.null(par) ){
    par <- c( 1, rep(1,ncol(x1)) )
  }
  sigma_f <- par[1]
  l <- par[2:length(par)]
  
  d1 <- nrow(x1)
  d2 <- nrow(x2)
  res <- matrix(-1, nrow=d1, ncol=d2)  
  if( d1 == d2 ){
    for( i in seq_len(d1) ){
      xi <- x1[i,]
      for( j in i:d2 ){
        res[i,j] <- sigma_f * exp( -0.5 * sum( l * (xi - x2[j,])^2) )  
      }
    }
    res <- as.matrix(forceSymmetric(res))
  } else{
    for( i in seq_len(d1) ){
      xi <- x1[i,] 
      for( j in seq_len(d2) ){
        res[i,j] <- sigma_f * exp( -0.5 * sum( l * (xi - x2[j,])^2) ) 
      }
    }  
  }  
  
  res
} )

optimizeHyperParameters <- function( x, y, numOpt=1, verbose=FALSE ){
  
  assign( "K", NULL, envir=.GlobalEnv )
  assign( "KInv", NULL, envir=.GlobalEnv )
  assign( "alpha", NULL, envir=.GlobalEnv )
  assign( "it", 1, envir=.GlobalEnv )
  
  fn <- calcLogLikelihood
  
  gr <- cmpfun( calcGradLogLikelihood )
  
  randomInit <- log( c( var(y), rep( 1/NCOL(x), NCOL(x) ), 0.1*var(y) ) )
  bestPar <- exp( optim( randomInit, fn=fn, gr=gr, control=list(fnscale=-1, maxit=150, trace=verbose, REPORT=1), method="BFGS",
                         x=x, y=y, eps=eps, verbose=verbose )$par )  
  
  bestPar
}

gpRegression <- cmpfun( function( u=NULL, y=NULL, lu=NULL, ly=NULL, x=NULL, par=NULL, kernel="SE", linearMean=FALSE, 
                                  optimization="ML", numOpt=1, optimizeOrder=FALSE, normalize=FALSE, verbose=TRUE, ... ){
  
  kernel <- toupper( kernel )
  if( kernel == "LIN" ){
    generateRandomParams <<- generateRandomParams_LIN
    calcKernelMatrix <<- calcKernelMatrix_LIN
    calcRegularizationPenalty <<- calcRegularizationPenalty_LIN
    calcGradLogLikelihood <<- calcGradLogLikelihood_LIN
  } else if( kernel == "SE" ){
    generateRandomParams <<- generateRandomParams_SE
    calcKernelMatrix <<- calcKernelMatrix_SE  
    calcRegularizationPenalty <<- calcRegularizationPenalty_SE
    calcGradLogLikelihood <<- calcGradLogLikelihood_SE
  } else if( kernel == "SELIN" ){
    generateRandomParams <<- generateRandomParams_SELIN
    calcKernelMatrix <<- calcKernelMatrix_SELIN
    calcRegularizationPenalty <<- calcRegularizationPenalty_SELIN
    calcGradLogLikelihood <<- calcGradLogLikelihood_SELIN
  } else if( kernel == "SELINMULT" ){
    generateRandomParams <<- generateRandomParams_SELINMULT
    calcKernelMatrix <<- calcKernelMatrix_SELINMULT
    calcRegularizationPenalty <<- calcRegularizationPenalty_SELINMULT
    calcGradLogLikelihood <<- calcGradLogLikelihood_SELINMULT
  } else if( kernel == "PER" ){
    generateRandomParams <<- generateRandomParams_PER
    calcKernelMatrix <<- calcKernelMatrix_PER
    calcRegularizationPenalty <<- calcRegularizationPenalty_PER
    calcGradLogLikelihood <<- calcGradLogLikelihood_PER
  } else if( kernel == "PERVAR" ){
    generateRandomParams <<- generateRandomParams_PERVAR
    calcKernelMatrix <<- calcKernelMatrix_PERVAR
    calcRegularizationPenalty <<- calcRegularizationPenalty_PERVAR
    calcGradLogLikelihood <<- calcGradLogLikelihood_PERVAR
  } else if( kernel == "MATERN" ){
    generateRandomParams <<- generateRandomParams_MATERN
    calcKernelMatrix <<- calcKernelMatrix_MATERN
    calcRegularizationPenalty <<- calcRegularizationPenalty_MATERN
    calcGradLogLikelihood <<- calcGradLogLikelihood_MATERN
  } else if( kernel == "RQ" ){
    generateRandomParams <<- generateRandomParams_RQ
    calcKernelMatrix <<- calcKernelMatrix_RQ
    calcRegularizationPenalty <<- calcRegularizationPenalty_RQ
    calcGradLogLikelihood <<- calcGradLogLikelihood_RQ
  } else if( kernel == "SM" ){
    generateRandomParams <<- generateRandomParams_SM
    calcKernelMatrix <<- calcKernelMatrix_SM
    calcRegularizationPenalty <<- calcRegularizationPenalty_SM
    calcGradLogLikelihood <<- calcGradLogLikelihood_SM
  } else if( kernel == "DSE" ){
    generateRandomParams <<- generateRandomParams_DSE
    calcKernelMatrix <<- calcKernelMatrix_DSE
    calcRegularizationPenalty <<- calcRegularizationPenalty_DSE
    calcGradLogLikelihood <<- calcGradLogLikelihood_DSE
  } else if( kernel == "ADD" ){
    generateRandomParams <<- generateRandomParams_ADD
    calcKernelMatrix <<- calcKernelMatrix_ADD
    calcRegularizationPenalty <<- calcRegularizationPenalty_ADD
    calcGradLogLikelihood <<- calcGradLogLikelihood_ADD
  } else if( kernel == "ADDLIN" ){
    generateRandomParams <<- generateRandomParams_ADDLIN
    calcKernelMatrix <<- calcKernelMatrix_ADDLIN
    calcRegularizationPenalty <<- calcRegularizationPenalty_ADDLIN
    calcGradLogLikelihood <<- calcGradLogLikelihood_ADDLIN
  }
  
  if( is.null(x) ){    
    if( optimizeOrder ){
      if( verbose ){ print("Optimzing orders...", quote=FALSE) }
      res <- gpOptimizeOrder( u=u, y=y, kernel=kernel )
      lu <- res$lu
      ly <- res$ly
      if( verbose ){ print(paste("Optimized orders: lu=", lu, ", ly=", ly ), quote=FALSE) }
    }
    
    if( verbose ){ print("Formatting estimation data...", quote=FALSE) }
    formatted <- formatData( u, y, lu, ly )
    x <- formatted$x 
    y <- formatted$y
  }
  
  lsPar <- NULL
  if( linearMean ){
    xLS <- cbind( 1, x )
    lsPar <- chol2inv( chol( t(xLS) %*% xLS ) ) %*% t(xLS) %*% y    
    y <- y - xLS %*% lsPar  
  }
  
  normX <- NULL
  normY <- NULL
  if( normalize ){
    normX <- normalizeData(x)
    x <- normX$data
    normY <- normalizeData(y)
    y <- normY$data
  }
  
  if( is.null(par) ){    
    if( verbose ){ print("Optimizing hyperparameters...", quote=FALSE) }
    optimization <- toupper( optimization )
    par <- optimizeHyperParameters( x, y, numOpt=numOpt, verbose=verbose )  
    if( verbose ){ print(par) }
  }
  
  sigmaNoise <- par[length(par)] + eps
  par <- par[1:(length(par)-1)]
  
  K <- calcKernelMatrix( x, x, par ) + sigmaNoise*diag(nrow(x))
  KInv <- chol2inv(jitchol(K))
  
  if( verbose ){ print("Estimation concluded.", quote=FALSE) }
  list( par=c(par,sigmaNoise), KInv=KInv, eps=eps, lsPar=lsPar,
        dataX=x, dataY=as.matrix(y), kernel=kernel, optimization=optimization, numOpt=numOpt,
        normX=normX, normY=normY, lu=lu, ly=ly )
} )

gpPredict <- cmpfun( function( m, u=NULL, y=NULL, x=NULL, fullPrediction=FALSE, verbose=TRUE, ... ){
  
  kernel <- toupper( m$kernel )
  if( kernel == "LIN" ){
    calcKernelMatrix <- calcKernelMatrix_LIN
  } else if( kernel == "SE" ){
    calcKernelMatrix <- calcKernelMatrix_SE  
  } else if( kernel == "SELIN" ){
    calcKernelMatrix <- calcKernelMatrix_SELIN
  } else if( kernel == "SELINMULT" ){
    calcKernelMatrix <- calcKernelMatrix_SELINMULT
  } else if( kernel == "PER" ){
    calcKernelMatrix <- calcKernelMatrix_PER
  } else if( kernel == "PERVAR" ){
    calcKernelMatrix <- calcKernelMatrix_PERVAR
  } else if( kernel == "MATERN" ){
    calcKernelMatrix <- calcKernelMatrix_MATERN
  } else if( kernel == "RQ" ){
    calcKernelMatrix <- calcKernelMatrix_RQ
  } else if( kernel == "SM" ){
    calcKernelMatrix <- calcKernelMatrix_SM
  } else if( kernel == "DSE" ){
    calcKernelMatrix <- calcKernelMatrix_DSE
  } else if( kernel == "ADD" ){
    calcKernelMatrix <- calcKernelMatrix_ADD
  } else if( kernel == "ADDLIN" ){
    calcKernelMatrix <- calcKernelMatrix_ADDLIN
  }
  
  lsPar <- m$lsPar
  
  if( verbose ){ print("Formatting test data...", quote=FALSE) }
  if( is.null(x) ){
    lu <- m$lu
    ly <- m$ly
    formatted <- formatData( u, y, lu, ly )
    x <- formatted$x
    y <- formatted$y
    if( ! is.null(m$normX) ){
      x <- normalizeData( x, mean=m$normX$mean, std=m$normX$std )$data
    }
    if( ! is.null(m$normY) ){
      y <- normalizeData( y, mean=m$normY$mean, std=m$normY$std )$data
    }
  } 
  
  sigmaNoise <- m$par[length(m$par)] 
  par <- m$par[1:(length(m$par)-1)]  
  KInv <- m$KInv
  dataX <- m$dataX
  dataY <- m$dataY
  
  if( verbose ){ print("Performing prediction...", quote=FALSE) }
  if( fullPrediction == FALSE ){
    yEst <- rep(0, nrow(x))
    yVar <- rep(0, nrow(x))
    auxMat <- KInv %*% dataY
    for( i in seq_len(nrow(x)) ){    
      KNewX <- calcKernelMatrix( x[i,,drop=FALSE], dataX, par )
      KNewNew <- calcKernelMatrix( x[i,,drop=FALSE], x[i,,drop=FALSE], par )
      
      yEst[i] <- as.numeric( KNewX %*% auxMat )
      yVar[i] <- as.numeric( KNewNew - KNewX %*% KInv %*% t(KNewX) )    
    }  
    yVar <- yVar + sigmaNoise
    
    if( is.null(lsPar) == FALSE ){
      yEst <- yEst + cbind( 1, x ) %*% lsPar  
    }
    
    if( ! is.null(m$normY) ){
      yEst <- yEst * m$normY$std + m$normY$mean
      yVar <- yVar * (m$normY$std^2)
    }
    
  } else{
    
    KNewX <- calcKernelMatrix( as.matrix(x), dataX, par )
    KNewNew <- calcKernelMatrix( as.matrix(x), as.matrix(x), par )
    
    yEst <- KNewX %*% KInv %*% dataY
    yVar <- KNewNew - KNewX %*% KInv %*% t(KNewX) + sigmaNoise*diag(NROW(x))
  }
  
  if( verbose ){ print("Prediction concluded.", quote=FALSE) }
  list( est=yEst, var=yVar )  
} )

gpSimulation <- cmpfun( function( m, u=NULL, y=NULL, verbose=TRUE, ... ){
  
  kernel <- toupper( m$kernel )
  if( kernel == "LIN" ){
    calcKernelMatrix <- calcKernelMatrix_LIN
  } else if( kernel == "SE" ){
    calcKernelMatrix <- calcKernelMatrix_SE  
  } else if( kernel == "SELIN" ){
    calcKernelMatrix <- calcKernelMatrix_SELIN
  } else if( kernel == "SELINMULT" ){
    calcKernelMatrix <- calcKernelMatrix_SELINMULT
  } else if( kernel == "PER" ){
    calcKernelMatrix <- calcKernelMatrix_PER
  } else if( kernel == "PERVAR" ){
    calcKernelMatrix <- calcKernelMatrix_PERVAR
  } else if( kernel == "MATERN" ){
    calcKernelMatrix <- calcKernelMatrix_MATERN
  } else if( kernel == "RQ" ){
    calcKernelMatrix <- calcKernelMatrix_RQ
  } else if( kernel == "SM" ){
    calcKernelMatrix <- calcKernelMatrix_SM
  } else if( kernel == "DSE" ){
    calcKernelMatrix <- calcKernelMatrix_DSE
  } else if( kernel == "ADD" ){
    calcKernelMatrix <- calcKernelMatrix_ADD
  } else if( kernel == "ADDLIN" ){
    calcKernelMatrix <- calcKernelMatrix_ADDLIN
  }
  
  if( verbose ){ print("Formatting test data...", quote=FALSE) }
  
  lu <- m$lu
  ly <- m$ly
  
  lsPar <- m$lsPar
  sigmaNoise <- m$par[length(m$par)] 
  par <- m$par[1:(length(m$par)-1)]  
  KInv <- m$KInv
  dataX <- m$dataX
  dataY <- m$dataY
  
  if( verbose ){ print("Performing free simulation...", quote=FALSE) }
  l <- max(lu,ly)
  yEst <- c( y[1:l], rep(0, length(u)-l) )
  yVar <- rep(sigmaNoise, length(u))
  uRegres <- rep(0, lu)
  yRegres <- rep(0, ly)
  
  xNew <- matrix( 0, nrow=length(u), ncol=lu+ly )
  
  auxMat <- KInv %*% dataY
  for( i in (l+1):length(u) ){
    
    iniU <- i-1; endU <- i-lu;
    iniY <- i-1; endY <- i-ly;
    uRegres <- u[iniU:endU] 
    yRegres <- yEst[iniY:endY]
    
    if( lu == 0 ){
      x <- t( as.matrix( c( yRegres ) ) )
    } else if( ly == 0 ){
      x <- t( as.matrix( uRegres ) )
    } else{
      x <- t( as.matrix( c( yRegres, uRegres ) ) )  
    }    
    
    if( ! is.null(m$normX) ){
      x <- normalizeData( x, mean=m$normX$mean, std=m$normX$std )$data
    }
    
    xNew[i,] <- x
    
    KNewX <- calcKernelMatrix( x, dataX, par )  
    KNewNew <- calcKernelMatrix( x, x, par )
    
    yEst[i] <- as.numeric( KNewX %*% auxMat )
    yVar[i] <- as.numeric( KNewNew - KNewX %*% KInv %*% t(KNewX) ) + sigmaNoise 
    
    if( is.null(lsPar) == FALSE ){
      yEst[i] <- yEst[i] + cbind( 1, x ) %*% lsPar  
    }
    
    if( ! is.null(m$normY) ){
      yEst[i] <- yEst[i] * m$normY$std + m$normY$mean
      yVar[i] <- yVar[i] * (m$normY$std^2)
    }
  }  
  
  if( verbose ){ print("Free simulation concluded.", quote=FALSE) }
  list( est=yEst, var=yVar, xNew=xNew )  
} )
