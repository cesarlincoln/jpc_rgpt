library(compiler)
library(Matrix)
source('code/gpFunctions.R', chdir=TRUE)
source('code/gpKernels.R', chdir=TRUE)

eps <- 10^-3

generateRandomParams <- NULL
calcGradLowerBound <- NULL

model <- NULL

calcPsi1 <- cmpfun( function( layer ){
  N <- NCOL(layer$muMat)
  M <- NCOL(layer$zMat)
  Q <- NROW(layer$zMat)
  psi1 <- matrix(0, nrow=N, ncol=M)
  for( n in 1:N ){    
    auxCalc1 <- layer$w * layer$sMat[,n] + 1
    if( any( is.na(auxCalc1) ) || any( auxCalc1 <= 0 ) ){ return(NA) }
    auxCalc2 <- 0.5*log(auxCalc1)
    psi1[n,] <- (layer$sigmaF/layer$sigma) * exp( .colSums( - 0.5*layer$w*(layer$muMat[,n] - layer$zMat)^2 / auxCalc1 - auxCalc2, Q, M ) )
  }
  psi1
} )
calcPsi2 <- cmpfun( function( layer ){
  N <- NCOL(layer$muMat)
  M <- NCOL(layer$zMat)
  Q <- NROW(layer$zMat)
  w <- layer$w
  coef <- (layer$sigmaF^2/layer$sigma) * exp( .colSums( - 0.25 * layer$w * layer$auxZMat1, Q, M*M ) )
  dPsi2ArrayCalc <- replicate( N, matrix(0, nrow=M, ncol=M), simplify=FALSE )
  for( n in 1:N ){    
    auxCalc1 <- 2 * layer$w * layer$sMat[,n] + 1
    auxCalc2 <- 0.5*log(auxCalc1)
    dPsi2ArrayCalc[[n]][,] <- coef *   
      exp( .colSums( ( - layer$w / auxCalc1) * (layer$muMat[,n] - 0.5*layer$auxZMat2)^2 - auxCalc2, Q, M*M ) )
  }
  dPsi2ArrayCalc
} )

calcLowerBound <- function( par, u, y, lu, ly, lag, parSkeleton, optimOpts=NULL, eps=0, verbose=FALSE ){
  
  yInit <- y[1:ly]
  y <- y[-(1:ly)]
  model <- relist(par, parSkeleton)
  N <- length(y)
  M <- NCOL(model[[1]]$zMat)
  numHiddenLayers <- length(model) - 1
  
  # Remove the logs
  
  for( h in 1:(numHiddenLayers+1) ){
    model[[h]]$zMat <- unname( model[[h]]$zMat )
    for( field in c( "sigma", "lambdaVarPar", "sigma0", "parF", "eps", "aVarPar", "bVarPar", "aPrior", "bPrior" ) ){
      if( any( is.null(model[[h]][[field]]) == FALSE ) ){
        model[[h]][[field]] <- exp( model[[h]][[field]] )
        if( field == "eps" ){
          model[[h]][[field]] <- pmax( model[[h]][[field]], 10^-10 )
        }
        if( field == "aVarPar" || field == "bVarPar" ){ 
          model[[h]][[field]] <- pmax( model[[h]][[field]], 10^-6 )
        }
      }
    }
  }
  
  # Create the variational parameters
  
  for( h in 1:(numHiddenLayers+1) ){
    layer <- model[[h]]
    muMat <- matrix( 0, nrow=NROW(layer$zMat), ncol=N )
    sMat <- matrix( 0, nrow=NROW(layer$zMat), ncol=N )
    
    if( h == numHiddenLayers+1 ){
      for( n in 1:N ){
        ind <- (n+ly-1):n
        muMat[,n] <- model[[h-1]]$muVarPar[1+ind]
        sMat[,n] <- model[[h-1]]$lambdaVarPar[1+ind] 
      }
      if( optimOpts[[h]]$robust ){
        model[[h]]$aVarPar <- model[[h]]$aVarPar 
        model[[h]]$bVarPar <- model[[h]]$bVarPar
        model[[h]]$R <- diag(model[[h]]$aVarPar/model[[h]]$bVarPar)
        model[[h]]$sigma <- 1
      }
    } else if( h == 1 ){
      for( n in 1:N ){
        ind <- (n+ly-1):n
        if( lu > 0 ){
          indu <- (ly-lu) + (n-lag+lu):(n-lag+1)
          muMat[,n] <- c( layer$muVarPar[ind], as.numeric(u[indu,]) )  
        } else{
          muMat[,n] <- layer$muVarPar[ind]
        }
        sMat[,n] <- c( layer$lambdaVarPar[ind], rep(0, lu*NCOL(u))) 
      }
    } else {
      for( n in 1:N ){
        ind <- (n+ly-1):n
        muMat[,n] <- c( layer$muVarPar[ind], model[[h-1]]$muVarPar[1+ind] )  
        sMat[,n] <- c( layer$lambdaVarPar[ind], model[[h-1]]$lambdaVarPar[1+ind] ) 
      }
    }
    
    model[[h]]$muMat <- muMat
    model[[h]]$sMat <- sMat
    model[[h]]$sigmaF <- layer$parF[1]
    model[[h]]$w <- layer$parF[-1]
  }
  
  # Calculate the covariance matrices
  
  for( h in 1:(numHiddenLayers+1) ){
    model[[h]]$Kmm <- calcKernelMatrix_SE( t(model[[h]]$zMat), t(model[[h]]$zMat), model[[h]]$parF ) + model[[h]]$eps*diag(M)
    model[[h]]$cholKmm <- jitchol( model[[h]]$Kmm, text="Kmm" )
    if( length(model[[h]]$cholKmm) == 1 && is.na(model[[h]]$cholKmm) ){ return(NA) }  
    model[[h]]$KmmInv <- chol2inv(model[[h]]$cholKmm)
  }
  
  # Calculate the PSI statistics
  
  for( h in 1:(numHiddenLayers+1) ){
    layer <- model[[h]]
    
    model[[h]]$psi0 <- N * layer$sigmaF
    
    model[[h]]$psi1 <- calcPsi1( model[[h]] )
    if( length(model[[h]]$psi1) == 1 && is.na(model[[h]]$psi1) ){ return(NA) }
    if( optimOpts[[h]]$robust ){
      model[[h]]$psi1 <- model[[h]]$R %*% model[[h]]$psi1  
    }
    
    model[[h]]$auxZMat1 <- (layer$zMat[,rep(1:M,M),drop=FALSE] - layer$zMat[,rep(1:M,each=M),drop=FALSE])^2
    model[[h]]$auxZMat2 <- layer$zMat[,rep(1:M,M),drop=FALSE] + layer$zMat[,rep(1:M,each=M),drop=FALSE]
    model[[h]]$auxZMat3 <- layer$zMat[,rep(1:M,M),drop=FALSE] - layer$zMat[,rep(1:M,each=M),drop=FALSE]
    model[[h]]$dPsi2ArrayCalc <- calcPsi2( model[[h]] )
    if( length(model[[h]]$psi2) == 1 && is.na(model[[h]]$psi2) ){ return(NA) }
    if( optimOpts[[h]]$robust ){
      model[[h]]$dPsi2ArrayCalc <- lapply( 1:N, function(n) model[[h]]$R[n,n] * model[[h]]$dPsi2ArrayCalc[[n]] )
      model[[h]]$psi2 <- Reduce( "+", model[[h]]$dPsi2ArrayCalc )
    } else{
      model[[h]]$psi2 <- Reduce( "+", model[[h]]$dPsi2ArrayCalc )  
    }
    
    At <- diag(M) + t( forwardsolve( t(layer$cholKmm), t(forwardsolve( t(layer$cholKmm), model[[h]]$psi2 )) ) )
    cholAt <- jitchol( At, text="At" )
    if( length(cholAt) == 1 && is.na(cholAt) ){ return(NA) }   
    model[[h]]$cholAt <- cholAt
    model[[h]]$auxMatInv <- t( backsolve( layer$cholKmm, t(backsolve( layer$cholKmm, chol2inv(cholAt) )) ) )
  }
  
  # Calculate the lower bound
  L <- 0
  for( h in 1:(numHiddenLayers+1) ){
    layer <- model[[h]]
    L <- L + 0.5*sum(diag( layer$KmmInv %*% layer$psi2 )) - sum(log(diag(layer$cholAt)))
    if( h < numHiddenLayers+1 ){
      L <- L - 0.5*N*log(2*pi*layer$sigma) +
        - (0.5/layer$sigma)*sum( layer$lambdaVarPar[(ly+1):(N+ly)] + layer$muVarPar[(ly+1):(N+ly)]^2 + layer$sigmaF ) +
        0.5*as.numeric( t(layer$muVarPar[(ly+1):(N+ly)]) %*% layer$psi1 %*% layer$auxMatInv %*% t(layer$psi1) %*% layer$muVarPar[(ly+1):(N+ly)] ) +
        0.5*(N+ly)*log(2*pi) + 0.5*(N+ly) + 0.5*sum( log(layer$lambdaVarPar) ) +
        - 0.5*ly*log(2*pi) - 0.5*sum(log(layer$sigma0)) - sum( (0.5/layer$sigma0) * ( layer$lambdaVarPar[1:ly] + layer$muVarPar[1:ly]^2 - 2*layer$muVarPar[1:ly]*layer$mu0 + layer$mu0^2 ) )
      
    } else{
      if( optimOpts[[h]]$robust ){
        L <- L - 0.5*N*log(2*pi) + 0.5*sum( digamma(layer$aVarPar) - log(layer$bVarPar) ) - 0.5*sum( diag(layer$R)*( y^2 + layer$sigmaF ) ) +
          - sum( (layer$aVarPar - layer$aPrior)*digamma(layer$aVarPar) - lgamma(layer$aVarPar) + lgamma(layer$aPrior) + layer$aPrior*(log(layer$bVarPar) - log(layer$bPrior)) + layer$aVarPar*(layer$bPrior - layer$bVarPar)/layer$bVarPar )  
      } else{
        L <- L - 0.5*N*log(2*pi*layer$sigma) - (0.5/layer$sigma)*sum( y^2 + layer$sigmaF )
      }
      L <- L + 0.5*as.numeric( t(y) %*% layer$psi1 %*% layer$auxMatInv %*% t(layer$psi1) %*% y )  
    }
  }
  
  model <<- model
  
  as.numeric(L)
}

calcGradLowerBound_SE <- function( par, u, y, lu, ly, lag, parSkeleton, optimOpts=NULL, eps=0, verbose=FALSE ){
  
  yInit <- y[1:ly]
  y <- y[-(1:ly)]
  
  if( ly > 1 ){
    par(mfrow=c(2,2))
    p1Z <- predict( prcomp( t( model[[1]]$zMat ) ), t( model[[1]]$zMat ) )[,1:2]
    p1X <- predict( prcomp( t( model[[1]]$muMat ) ), t( model[[1]]$muMat ) )[,1:2]
    plot( p1Z, col="red", pch="*", cex=2.5,
          xlim=range( c( p1Z[,1], p1X[,1] ) ),
          ylim=range( c( p1Z[,2], p1X[,2] ) ) )
    points( p1X )
    p2Z <- predict( prcomp( t( model[[2]]$zMat ) ), t( model[[2]]$zMat ) )[,1:2]
    p2X <- predict( prcomp( t( model[[2]]$muMat ) ), t( model[[2]]$muMat ) )[,1:2]
    plot( p2Z, col="red", pch="*", cex=2.5,
          xlim=range( c( p2Z[,1], p2X[,1] ) ),
          ylim=range( c( p2Z[,2], p2X[,2] ) ) )
    points( p2X )
    if( numHiddenLayers > 1 ){
      p3Z <- predict( prcomp( t( model[[3]]$zMat ) ), t( model[[3]]$zMat ) )[,1:2]
      p3X <- predict( prcomp( t( model[[3]]$muMat ) ), t( model[[3]]$muMat ) )[,1:2]
      plot( p3Z, col="red", pch="*", cex=2.5,
            xlim=range( c( p3Z[,1], p3X[,1] ) ),
            ylim=range( c( p3Z[,2], p3X[,2] ) ) )
      points( p3X )
    }
  }
  plot( c(yInit, y), type="l", ylab="" )
  colorVec <- c( palette()[2:(numHiddenLayers+1)], "black" )
  for( h in 1:(numHiddenLayers+1) ){
    # cat("Layer", h, "\n")
    # cat("parF", model[[h]]$parF, "\n")
    # cat("sigma", model[[h]]$sigma, "\n")
    if( h <= numHiddenLayers ){
      lines( model[[h]]$muVarPar, col=colorVec[h] )
      # cat("muVarPar", range(model[[h]]$muVarPar), "\n")
      # cat("lambdaVarPar", range(model[[h]]$lambdaVarPar), "\n")
    }
    # cat("\n")
  }
  legend( "topleft", legend=c( sapply( 1:numHiddenLayers, function(h) paste("Layer", h) ), "Training output" ), col=colorVec, lwd=2 )
  
  # Outliers location
  # plot( diag(model[[numHiddenLayers+1]]$R) )
  # points( d$outlierIndex-1, diag(model[[numHiddenLayers+1]]$R)[d$outlierIndex-1], col="red", pch="x" )
  
  grPar <- relist(par, parSkeleton)
  N <- length(y)
  M <- NCOL(model[[1]]$zMat)
  numHiddenLayers <- length(model) - 1
  yyt <- tcrossprod(y)
  
  for( h in 1:(numHiddenLayers+1) ){
    
    # Pre-computations
    
    layer <- model[[h]]
    Q <- length(layer$parF[-1])
    zeroMatMM <- matrix( 0, nrow=M, ncol=M )
    zeroMatNM <- matrix( 0, nrow=N, ncol=M )
    
    # Gradients wrt kernel hyperparameters
    
    if( optimOpts[[h]]$holdEps == FALSE ){
      dKmmEps <- diag(M)
      
      grEps <- 0.5*sum(diag( - layer$KmmInv %*% dKmmEps %*% layer$KmmInv %*% layer$psi2 )) +
        0.5*sum(diag( layer$KmmInv %*% dKmmEps )) - 0.5*sum(diag( layer$auxMatInv %*% ( dKmmEps ) ))
      
      if( h < numHiddenLayers+1 ){
        grEps <- grEps + 0.5 * t(layer$muVarPar[(ly+1):(N+ly)]) %*% ( - layer$psi1 %*% layer$auxMatInv %*% dKmmEps %*% layer$auxMatInv %*% t(layer$psi1) ) %*% layer$muVarPar[(ly+1):(N+ly)] 
      } else{
        grEps <- grEps + 0.5 * t(y) %*% ( - layer$psi1 %*% layer$auxMatInv %*% dKmmEps %*% layer$auxMatInv %*% t(layer$psi1) ) %*% y
      }
      grEps <- grEps * layer$eps
    } else{
      grEps <- 0
    }
    
    if( optimOpts[[h]]$holdParF == FALSE ){
      if( optimOpts[[h]]$holdFirst == FALSE ){
        dPsi1SigmaF <- layer$psi1 / layer$sigmaF
        dPsi2SigmaF <- 2 * layer$psi2 / layer$sigmaF
        dKmmSigmaF <- ( layer$Kmm - layer$eps*diag(M) ) / layer$sigmaF
        
        if( h == numHiddenLayers+1 && optimOpts[[h]]$robust ){
          grSigmaF <- - 0.5*sum(diag(layer$R))
        } else{
          grSigmaF <- - 0.5*N/layer$sigma  
        }
        
        grSigmaF <- grSigmaF + 0.5*sum(diag( - layer$KmmInv %*% dKmmSigmaF %*% layer$KmmInv %*% layer$psi2 + layer$KmmInv %*% dPsi2SigmaF )) +
          0.5*sum(diag( layer$KmmInv %*% dKmmSigmaF )) - 0.5*sum(diag( layer$auxMatInv %*% ( dKmmSigmaF + dPsi2SigmaF ) ))
        
        if( h < numHiddenLayers+1 ){
          grSigmaF <- grSigmaF + 0.5 * t(layer$muVarPar[(ly+1):(N+ly)]) %*% ( dPsi1SigmaF %*% layer$auxMatInv %*% t(layer$psi1) - layer$psi1 %*% layer$auxMatInv %*% ( dKmmSigmaF + dPsi2SigmaF ) %*% layer$auxMatInv %*% t(layer$psi1) +
                                                                                layer$psi1 %*% layer$auxMatInv %*% t(dPsi1SigmaF) ) %*% layer$muVarPar[(ly+1):(N+ly)] 
        } else{
          grSigmaF <- grSigmaF + 0.5 * t(y) %*% ( dPsi1SigmaF %*% layer$auxMatInv %*% t(layer$psi1) - layer$psi1 %*% layer$auxMatInv %*% ( dKmmSigmaF + dPsi2SigmaF ) %*% layer$auxMatInv %*% t(layer$psi1) +
                                                    layer$psi1 %*% layer$auxMatInv %*% t(dPsi1SigmaF) ) %*% y
        }
        grSigmaF <- grSigmaF * layer$sigmaF
      } else{
        grSigmaF <- 0  
      }
      
      grW <- rep(0, Q)
      aux1 <- tcrossprod( layer$auxMatInv, layer$psi1 )
      if( h < numHiddenLayers+1 ){
        aux2 <- tcrossprod( layer$muVarPar[(ly+1):(N+ly)] )
        aux3 <- aux1 %*% aux2 %*% t(aux1)
        auxCalcW1 <- 2 * t( aux1 %*% aux2 )
      } else{
        aux3 <- aux1 %*% yyt %*% t(aux1)
        auxCalcW1 <- 2 * t( aux1 %*% yyt )
      }
      auxCalcW2 <- - ( aux3 + layer$auxMatInv + layer$KmmInv %*% layer$psi2 %*% layer$KmmInv - layer$KmmInv )
      auxCalcW3 <- layer$KmmInv - layer$auxMatInv - aux3
      for( q in 1:Q ){
        
        dPsi1WMat <- matrix( 0, nrow=N, ncol=M )
        for( m in 1:M ){
          dPsi1WMat[,m] <- layer$psi1[,m] * (
            ( - 0.5 * (layer$muMat[q,] - layer$zMat[q,m] )^2 ) * ( 1/(layer$w[q]*layer$sMat[q,] + 1) - layer$w[q] * layer$sMat[q,]/(layer$w[q]*layer$sMat[q,] + 1)^2 ) +
              - 0.5*layer$sMat[q,]/(layer$w[q]*layer$sMat[q,] + 1) )
        }
        
        dPsi2WMat <- ( - 0.25 * layer$auxZMat1[q,] ) * layer$psi2
        auxCalc1 <- 1/(2*layer$w[q]*layer$sMat[q,] + 1) - layer$w[q] * 2*layer$sMat[q,]/(2*layer$w[q]*layer$sMat[q,] + 1)^2
        auxCalc2 <- layer$sMat[q,]/(2*layer$w[q]*layer$sMat[q,] + 1)
        # for( n in 1:N ){
        #   dPsi2WMat <- dPsi2WMat + layer$dPsi2ArrayCalc[[n]] *
        #     ( (layer$muMat[q,n] - 0.5*layer$auxZMat2[q,])^2 * ( - auxCalc1[n] ) - auxCalc2[n] )
        # }
        aux1 <- 0.5*layer$auxZMat2[q,]
        for( n in 1:N ){
          dPsi2WMat <- dPsi2WMat + layer$dPsi2ArrayCalc[[n]] *
            ( (layer$muMat[q,n] - aux1)^2 * ( - auxCalc1[n] ) - auxCalc2[n] )
        }
        
        dKmmW <- - 0.5 * outer( layer$zMat[q,], layer$zMat[q,], function(p,q) (p-q)^2 ) * ( layer$Kmm - layer$eps*diag(M) )
        
        grW[q] <- 0.5 * sum( auxCalcW2 * dKmmW ) + 0.5 * sum( auxCalcW3 * dPsi2WMat ) + 0.5 * sum( auxCalcW1 * dPsi1WMat )
      }
      grW <- grW * layer$w
      
      grParF <- c( grSigmaF, grW )
      
    } else{
      grParF <- rep( 0, 1+Q )
    }
    
    # Gradients wrt noise variance
    
    if( optimOpts[[h]]$holdSigma == FALSE ){
      dPsi1Sigma <- - layer$psi1 / layer$sigma
      dPsi2Sigma <- - layer$psi2 / layer$sigma
      grSigma <- - 0.5*N/layer$sigma +
        0.5*sum(diag( layer$KmmInv %*% dPsi2Sigma )) - 0.5*sum(diag( layer$auxMatInv %*% dPsi2Sigma ))
      
      if( h < numHiddenLayers+1 ){
        grSigma <- grSigma + (0.5/layer$sigma^2)*sum( layer$lambdaVarPar[(ly+1):(N+ly)] + layer$muVarPar[(ly+1):(N+ly)]^2 + layer$sigmaF ) +
          0.5*t(layer$muVarPar[(ly+1):(N+ly)]) %*% ( 2 * dPsi1Sigma %*% layer$auxMatInv %*% t(layer$psi1) - layer$psi1 %*% layer$auxMatInv %*% dPsi2Sigma %*% layer$auxMatInv %*% t(layer$psi1) ) %*% layer$muVarPar[(ly+1):(N+ly)]
      } else{
        if( optimOpts[[h]]$robust ){
          grSigma <- 0  
        } else{
          grSigma <- grSigma + (0.5/layer$sigma^2)*sum( y*y + layer$sigmaF ) +
            0.5*t(y) %*% ( 2 * dPsi1Sigma %*% layer$auxMatInv %*% t(layer$psi1) - layer$psi1 %*% layer$auxMatInv %*% dPsi2Sigma %*% layer$auxMatInv %*% t(layer$psi1) ) %*% y  
        }
      }
      
      grSigma <- grSigma * layer$sigma
    } else{
      grSigma <- 0
    }
    
    # Gradients wrt pseudo-inputs
    
    if( optimOpts[[h]]$holdZ == FALSE ){
      
      # dPsi1ZArray <- replicate( M, matrix(0, nrow=N, ncol=Q), simplify=FALSE )
      # for( n in 1:N ){
      #   dPsi1ZAux <- t( layer$psi1[n,] * t( layer$w * ( layer$muMat[,n] - layer$zMat ) / (layer$w*layer$sMat[,n] + 1) ) )
      #   for( m in 1:M ){
      #     dPsi1ZArray[[m]][n,] <- dPsi1ZAux[,m]
      #   }
      # }
      auxList <- lapply( 1:N, function(n) layer$psi1[n,] * t( layer$w * ( layer$muMat[,n] - layer$zMat ) / (layer$w*layer$sMat[,n] + 1) ) )
      dPsi1ZArray <- lapply( 1:M, function(m) matrix( t( sapply( auxList, function(l) l[m,] ) ), N, Q ) )
      
      # dPsi2ZArray <- lapply( 1:Q, function(q) - layer$w[q] * layer$auxZMat3[q,] * layer$psi2 )
      # auxCalc <- 2 * layer$w/(2*layer$w*layer$sMat + 1)
      # for( q in 1:Q ){
      #   for( n in 1:N ){
      #     dPsi2ZArray[[q]] <- dPsi2ZArray[[q]] + layer$dPsi2ArrayCalc[[n]] * auxCalc[q,n] * ( layer$muMat[q,n] - 0.5*layer$auxZMat2[q,] )  
      #   }
      # }
      dPsi2ZArray <- lapply( 1:Q, function(q) - layer$w[q] * layer$auxZMat3[q,] * layer$psi2 )
      auxCalc <- 2 * layer$w/(2*layer$w*layer$sMat + 1)
      for( q in 1:Q ){
        aux1 <- auxCalc[q,]
        aux2 <- layer$muMat[q,]
        aux3 <- - 0.5*layer$auxZMat2[q,]
        for( n in 1:N ){
          dPsi2ZArray[[q]] <- dPsi2ZArray[[q]] + layer$dPsi2ArrayCalc[[n]] * aux1[n] * ( aux2[n] + aux3 )
        }
      }
      
      grZMat <- matrix( 0, nrow=Q, ncol=M )
      aux1 <- layer$psi1 %*% layer$auxMatInv 
      if( h < numHiddenLayers + 1 ){
        aux2 <- tcrossprod( layer$muVarPar[(ly+1):(N+ly)] )
        aux3 <- crossprod( t(layer$muVarPar[(ly+1):(N+ly)]) %*% aux1 )  
      } else{
        aux2 <- yyt
        aux3 <- crossprod( t(y) %*% aux1 )  
      }
      auxcalcZ <- layer$KmmInv - layer$auxMatInv - layer$KmmInv %*% layer$psi2 %*% layer$KmmInv - aux3
      auxcalcZ2 <- layer$auxMatInv - layer$KmmInv + aux3
      auxcalcZ3 <- 2*aux2 %*% aux1
      auxKmm <- layer$Kmm - layer$eps*diag(M)
      for( m in 1:M ){
        for( q in 1:Q ){
          dKmmMat <- dPsi2ZMat <- zeroMatMM
          dKmmMat[m,] <- - auxKmm[m,] * layer$w[q] * ( layer$zMat[q,m] - layer$zMat[q,] )
          dKmmMat[,m] <- dKmmMat[m,]
          dPsi1ZMat <- zeroMatNM
          dPsi1ZMat[,m] <- dPsi1ZArray[[m]][,q]
          dPsi2ZMat[m,] <- dPsi2ZArray[[q]][m,]
          grZMat[q,m] <- 0.5*sum( auxcalcZ * dKmmMat ) - 0.5*sum( auxcalcZ2 * dPsi2ZMat ) +
            0.5 * sum( auxcalcZ3 * dPsi1ZMat )
        }
      }
      grZMat <- as.numeric( grZMat )
      
    } else{
      
      grZMat <- rep( 0, length(layer$zMat ) )  
    }
    
    if( h < numHiddenLayers+1 ){
      
      # Gradients wrt the variational parameters
      
      if( optimOpts[[h]]$holdX == FALSE ){
        
        if( h == numHiddenLayers ){
          indexNext <- 0
        } else{
          indexNext <- ly  
        }
        
        grMuVarPar <- rep(0,N+ly)
        aux1 <- layer$psi1 %*% layer$auxMatInv 
        aux2 <- model[[h+1]]$psi1 %*% model[[h+1]]$auxMatInv 
        auxCalcMu1 <- layer$KmmInv - layer$auxMatInv - crossprod( t(layer$muVarPar[(ly+1):(N+ly)]) %*% aux1 )
        auxCalcMu2 <- 2*tcrossprod(layer$muVarPar[(ly+1):(N+ly)]) %*% aux1
        if( h+1 <= numHiddenLayers ){
          auxCalcMu3 <- model[[h+1]]$KmmInv - model[[h+1]]$auxMatInv - crossprod( t(model[[h+1]]$muVarPar[(ly+1):(N+ly)]) %*% aux2 )
          auxCalcMu4 <- 2*tcrossprod(model[[h+1]]$muVarPar[(ly+1):(N+ly)]) %*% aux2  
        } else{
          auxCalcMu3 <- model[[h+1]]$KmmInv - model[[h+1]]$auxMatInv - crossprod( t(y) %*% aux2 )
          auxCalcMu4 <- 2*tcrossprod(y) %*% aux2
        }
        for( n in 1:(N+ly) ){
          
          dPsi1MuVarPar <- zeroMatNM
          dPsi2MuVarPar<- zeroMatMM
          for( q in (1:ly)[n+1:ly-ly > 0 & n+1:ly-ly <= N] ){
            dPsi1MuVarPar[n+q-ly,] <- dPsi1MuVarPar[n+q-ly,] +
              layer$psi1[n+q-ly,] * layer$w[q] * ( - (layer$muMat[q,n+q-ly] - layer$zMat[q,]) ) / (layer$w[q]*layer$sMat[q,n+q-ly] + 1)
            
            dPsi2MuVarPar <- dPsi2MuVarPar +
              layer$dPsi2ArrayCalc[[n+q-ly]] * (layer$muMat[q,n+q-ly] - 0.5*layer$auxZMat2[q,]) * ( - layer$w[q] * 2 / (2*layer$w[q]*layer$sMat[q,n+q-ly] + 1) )
          }
          
          dPsi1NextMuVarPar <- zeroMatNM
          dPsi2NextMuVarPar <- zeroMatMM
          for( q in (1:ly)[n+1:ly-ly-1 > 0 & n+1:ly-ly-1 <= N] ){
            dPsi1NextMuVarPar[n+q-ly-1,] <- dPsi1NextMuVarPar[n+q-ly-1,] +
              model[[h+1]]$psi1[n+q-ly-1,] * model[[h+1]]$w[q+indexNext] * ( - (model[[h+1]]$muMat[q+indexNext,n+q-ly-1] - model[[h+1]]$zMat[q+indexNext,]) ) / (model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1)  
            
            dPsi2NextMuVarPar <- dPsi2NextMuVarPar +
              model[[h+1]]$dPsi2ArrayCalc[[n+q-ly-1]] * (model[[h+1]]$muMat[q+indexNext,n+q-ly-1] - 0.5*model[[h+1]]$auxZMat2[q+indexNext,]) * ( - model[[h+1]]$w[q+indexNext] * 2 / (2*model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1) )  
          }
          
          grMuVarPar[n] <- 0.5*sum( auxCalcMu3 * dPsi2NextMuVarPar ) + 0.5*sum( auxCalcMu4 * dPsi1NextMuVarPar ) + 
            0.5*sum( auxCalcMu1 * dPsi2MuVarPar ) + 0.5*sum( auxCalcMu2 * dPsi1MuVarPar )
        }
        
        grMuVarPar[(ly+1):(N+ly)] <- grMuVarPar[(ly+1):(N+ly)] +
          - (1/layer$sigma)*layer$muVarPar[(ly+1):(N+ly)] + layer$psi1 %*% layer$auxMatInv %*% t(layer$psi1) %*% layer$muVarPar[(ly+1):(N+ly)]
        grMuVarPar[1:ly] <- grMuVarPar[1:ly] - ( (0.5/layer$sigma0) * ( 2*layer$muVarPar[1:ly] - 2*layer$mu0 ) )
        
        grLambdaVarPar <- rep(0,N+ly)
        for( n in 1:(N+ly) ){
          
          dPsi1LambdaVarPar <- zeroMatNM
          dPsi2LambdaVarPar <- zeroMatMM
          for( q in (1:ly)[n+1:ly-ly > 0 & n+1:ly-ly <= N] ){
            dPsi1LambdaVarPar[n+q-ly,] <- dPsi1LambdaVarPar[n+q-ly,] +
              layer$psi1[n+q-ly,] * ( 0.5*layer$w[q]^2 * (layer$muMat[q,n+q-ly] - layer$zMat[q,])^2 / (layer$w[q]*layer$sMat[q,n+q-ly] + 1)^2 +
                                        - 0.5*layer$w[q]/(layer$w[q]*layer$sMat[q,n+q-ly] + 1) )
            
            dPsi2LambdaVarPar <- dPsi2LambdaVarPar +
              layer$dPsi2ArrayCalc[[n+q-ly]] * ( (layer$muMat[q,n+q-ly] - 0.5*layer$auxZMat2[q,])^2 * ( 2*layer$w[q]^2 / (2*layer$w[q]*layer$sMat[q,n+q-ly] + 1)^2 ) +
                                                   - layer$w[q]/(2*layer$w[q]*layer$sMat[q,n+q-ly] + 1) )
          }
          
          dPsi1NextLambdaVarPar <- zeroMatNM
          dPsi2NextLambdaVarPar <- zeroMatMM
          for( q in (1:ly)[n+1:ly-ly-1 > 0 & n+1:ly-ly-1 <= N] ){
            dPsi1NextLambdaVarPar[n+q-ly-1,] <- dPsi1NextLambdaVarPar[n+q-ly-1,] +
              model[[h+1]]$psi1[n+q-ly-1,] * ( 0.5*model[[h+1]]$w[q+indexNext]^2 * (model[[h+1]]$muMat[q+indexNext,n+q-ly-1] - model[[h+1]]$zMat[q+indexNext,])^2 / (model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1)^2 +
                                                 - 0.5*model[[h+1]]$w[q+indexNext]/(model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1) ) 
            
            dPsi2NextLambdaVarPar <- dPsi2NextLambdaVarPar +
              model[[h+1]]$dPsi2ArrayCalc[[n+q-ly-1]] * ( (model[[h+1]]$muMat[q+indexNext,n+q-ly-1] - 0.5*model[[h+1]]$auxZMat2[q+indexNext,])^2 * ( 2*model[[h+1]]$w[q+indexNext]^2 / (2*model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1)^2 ) +
                                                            - model[[h+1]]$w[q+indexNext]/(2*model[[h+1]]$w[q+indexNext]*model[[h+1]]$sMat[q+indexNext,n+q-ly-1] + 1) )    
          }
          
          grLambdaVarPar[n] <- 0.5*sum( auxCalcMu3 * dPsi2NextLambdaVarPar ) + 0.5*sum( auxCalcMu4 * dPsi1NextLambdaVarPar ) + 
            0.5*sum( auxCalcMu1 * dPsi2LambdaVarPar ) + 0.5*sum( auxCalcMu2 * dPsi1LambdaVarPar )
        }
        
        grLambdaVarPar <- grLambdaVarPar + 0.5/layer$lambdaVarPar
        grLambdaVarPar[(ly+1):(N+ly)] <- grLambdaVarPar[(ly+1):(N+ly)] - (0.5/layer$sigma)
        grLambdaVarPar[1:ly] <- grLambdaVarPar[1:ly] - (0.5/layer$sigma0)
        grLambdaVarPar <- grLambdaVarPar * layer$lambdaVarPar
        
      } else{
        
        grMuVarPar <- rep( 0, length(layer$muVarPar) )
        grLambdaVarPar <- rep( 0, length(layer$lambdaVarPar) )
      }
      
      # Gradients wrt prior hyperparameters
      
      if( optimOpts[[h]]$holdPrior == FALSE ){
        grMu0 <- - (0.5/layer$sigma0) * ( - 2*layer$muVarPar[1:ly] + 2*layer$mu0 )
        
        grSigma0 <- - 0.5/layer$sigma0 + (0.5/layer$sigma0^2) * ( layer$lambdaVarPar[1:ly] + layer$muVarPar[1:ly]^2 - 2*layer$muVarPar[1:ly]*layer$mu0 + layer$mu0^2 )
        grSigma0 <- grSigma0 * layer$sigma0
        
      } else{
        grMu0 <- rep( 0, ly )  
        grSigma0 <- rep( 0, ly )  
      }
    } else{
      
      if( optimOpts[[h]]$robust ){
        
        # Gradients wrt variational parameters of the Student-t likelihood
        
        auxAB <- 0.5 * crossprod( t(y) %*% layer$psi1 %*% layer$auxMatInv )
        aux1 <- layer$auxMatInv %*% t(layer$psi1) %*% yyt
        if( optimOpts[[h]]$holdA == FALSE ){
          
          grAVarPar <- 0.5*trigamma(layer$aVarPar) - 0.5/layer$bVarPar*( y*y + layer$sigmaF ) +
            ( layer$aPrior - layer$aVarPar ) * trigamma( layer$aVarPar ) - ( layer$bPrior - layer$bVarPar ) / layer$bVarPar
          for( n in 1:N ){
            dPsi1AVarPar <- matrix( 0, nrow=M, ncol=N )
            dPsi1AVarPar[,n] <- layer$psi1[n,] / layer$aVarPar[n]
            aux2 <- sum( aux1 * dPsi1AVarPar )
            dPsi2AVarPar <- layer$dPsi2ArrayCalc[[n]] / layer$aVarPar[n]
            grAVarPar[n] <- grAVarPar[n] + 0.5*sum( layer$KmmInv * dPsi2AVarPar ) - 0.5*sum( layer$auxMatInv *  dPsi2AVarPar ) +
              aux2 - sum( auxAB * dPsi2AVarPar)
          }
          grAVarPar <- grAVarPar * layer$aVarPar
        } else{
          grAVarPar <- rep( 0, length(layer$aVarPar) )
        }
        
        if( optimOpts[[h]]$holdB == FALSE ){
          grBVarPar <- - 0.5/layer$bVarPar + (0.5*layer$aVarPar/layer$bVarPar^2)*( y*y + layer$sigmaF ) +
            ( - layer$aPrior / layer$bVarPar + layer$aVarPar * layer$bPrior / layer$bVarPar^2 ) 
          for( n in 1:N ){
            dPsi1BVarPar <- matrix( 0, nrow=M, ncol=N )
            dPsi1BVarPar[,n] <- - layer$psi1[n,] / layer$bVarPar[n]
            aux2 <- sum( aux1 * dPsi1BVarPar )
            dPsi2BVarPar <- - layer$dPsi2ArrayCalc[[n]] / layer$bVarPar[n]
            grBVarPar[n] <- grBVarPar[n] + 0.5*sum( layer$KmmInv * dPsi2BVarPar ) - 0.5*sum( layer$auxMatInv *  dPsi2BVarPar ) +
              aux2 - sum( auxAB * dPsi2BVarPar)
          }
          grBVarPar <- grBVarPar * layer$bVarPar
        } else{
          grBVarPar <- rep( 0, length(layer$bVarPar) )  
        }
        
        # Gradients wrt prior hyperparameters
        
        if( optimOpts[[h]]$holdPrior == FALSE ){
          if( optimOpts[[h]]$holdA == FALSE ){
            grAPrior <- - sum( - digamma(layer$aVarPar) + digamma(layer$aPrior) + log(layer$bVarPar) - log(layer$bPrior) )
            grAPrior <- grAPrior * layer$aPrior  
          } else{
            grAPrior <- 0  
          }
          
          if( optimOpts[[h]]$holdB == FALSE ){
            grBPrior <- - sum( - layer$aPrior/layer$bPrior + layer$aVarPar/layer$bVarPar )
            grBPrior <- grBPrior * layer$bPrior
          } else{
            grBPrior <- 0  
          }
          
        } else{
          grAPrior <- 0  
          grBPrior <- 0  
        }
      }
    }
    
    grPar[[h]]$parF <- grParF
    grPar[[h]]$sigma <- grSigma
    grPar[[h]]$eps <- grEps
    grPar[[h]]$zMat <- grZMat
    if( h < numHiddenLayers+1 ){
      grPar[[h]]$muVarPar <- grMuVarPar
      grPar[[h]]$lambdaVarPar <- grLambdaVarPar
      grPar[[h]]$mu0 <- grMu0
      grPar[[h]]$sigma0 <- grSigma0  
    } else{
      if( optimOpts[[h]]$robust ){
        grPar[[h]]$aVarPar <- grAVarPar
        grPar[[h]]$bVarPar <- grBVarPar
        grPar[[h]]$aPrior <- grAPrior
        grPar[[h]]$bPrior <- grBPrior  
      }
    }
  }
  
  as.numeric(unlist(grPar))
}

optimizeHyperParameters <- function( u, y, lu, ly, lag=1, M=10, numHiddenLayers=1, numParDynKernel=1, kernel=kernel, robustFlag=TRUE, verbose=FALSE ){
  
  fn <- cmpfun( calcLowerBound )
  gr <- cmpfun( calcGradLowerBound_SE )
  
  yInit <- y[1:ly]
  y <- y[-(1:ly)]
  if( lu > 0 ){
    uInit <- u[1:lu,,drop=FALSE]
    u <- u[-(1:lu),,drop=FALSE]
  }
  uTotal <- NULL
  if( lu > 0 ){
    uTotal <- unname( rbind( uInit, u ) )
  }
  
  # Initialization parameters
  N <- length(y)
  sigmaInit <- 0.01
  if( numHiddenLayers == 1 ){
    lambdaInit <- 0.5
  } else{
    lambdaInit <- 0.1
  }
  
  optimInit <- vector( "list", numHiddenLayers+1 )
  optimOpts <- vector( "list", numHiddenLayers+1 )
  library(cluster)
  for( h in 1:(numHiddenLayers+1) ){
    if( h == numHiddenLayers+1 ){
      xH <- matrix( 0, nrow=N, ncol=ly )
      for( n in 1:N ){
        ind <- (n+ly-1):n
        xH[n,] <- optimInit[[h-1]]$muVarPar[1+ind]
      }
      
      if( robustFlag ){
        optimInit[[h]]$aVarPar <- log( rep(1,N) ) 
        optimInit[[h]]$bVarPar <- log( rep(0.01,N) )
        optimInit[[h]]$aPrior <- log(1)
        optimInit[[h]]$bPrior <- log(0.01)
      } else{
        optimOpts[[h]]$robust <- FALSE
      }
      
    } else{
      dataInit <- as.matrix( c( yInit, y ) )
      if( lu > 0 ){
        dataInit <- cbind( dataInit, c( uInit, u ) ) 
      }
      dataInit <- unname( scale( dataInit ) )
      if( lu > 0 ){
        optimInit[[h]]$muVarPar <- unname( scale( predict( prcomp( dataInit ), dataInit )[,1] ) )
      } else{
        optimInit[[h]]$muVarPar <- unname( scale( c( yInit, y ) ) )
      }
      optimInit[[h]]$lambdaVarPar <- rep( log( lambdaInit * var( optimInit[[h]]$muVarPar ) ), N+ly )
      optimInit[[h]]$mu0 <- rep( 0, ly )
      optimInit[[h]]$sigma0 <- rep( log(1), ly )
      
      if( h == 1 ){
        xH <- matrix( 0, nrow=N, ncol=ly+lu*NCOL(u) )
        for( n in 1:N ){
          ind <- (n+ly-1):n
          if( lu > 0 ){
            indu <- (ly-lu) + (n-lag+lu):(n-lag+1)
            xH[n,] <- c( optimInit[[h]]$muVarPar[ind], as.numeric(uTotal[indu,]) )  
          } else{
            xH[n,] <- optimInit[[h]]$muVarPar[ind]
          }
        }
      } else {
        xH <- matrix( 0, nrow=N, ncol=2*ly )
        for( n in 1:N ){
          ind <- (n+ly-1):n
          xH[n,] <- c( optimInit[[h]]$muVarPar[ind], optimInit[[h-1]]$muVarPar[1+ind] )  
        }
      }  
    }
    
    optimInit[[h]]$sigma <- log( sigmaInit * ifelse( is.null(optimInit[[h]]$muVarPar), var(y), var(optimInit[[h]]$muVarPar) ) )   
    
    optimInit[[h]]$zMat <- t( pam(xH, M)$medoids )
    optimInit[[h]]$parF <- log( c( 1, rep( 1/NCOL(xH), NCOL(xH) ) ) )
    
    optimInit[[h]]$eps <- log( eps * exp(optimInit[[h]]$parF[1]) )
    
    if( h == numHiddenLayers+1 && robustFlag ){
      optimOpts[[h]]$robust <- TRUE  
    } else{
      optimOpts[[h]]$robust <- FALSE
    }
    optimOpts[[h]]$holdFirst <- FALSE
    optimOpts[[h]]$holdSigma <- FALSE
    optimOpts[[h]]$holdParF <- FALSE
    optimOpts[[h]]$holdEps <- FALSE
    optimOpts[[h]]$holdX <- FALSE
    optimOpts[[h]]$holdZ <- FALSE
    optimOpts[[h]]$holdPrior <- FALSE
    optimOpts[[h]]$holdA <- FALSE
    optimOpts[[h]]$holdB <- FALSE
  }
  names(optimInit) <- sapply( 1:(numHiddenLayers+1), function(h) paste("layer", h, sep="") )
  names(optimOpts) <- sapply( 1:(numHiddenLayers+1), function(h) paste("layer", h, sep="") )
  
  print( "Starting optimization...", quote=FALSE )
  
  # Settings for the optimization
  optimInit <- as.relistable( optimInit )
  
  randomInit <- unlist( optimInit )
  if( robustFlag ){
    
    optimOpts <- lapply( optimOpts, "[[<-", "holdSigma", TRUE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdFirst", TRUE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdEps", TRUE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdPrior", TRUE )
    optimOpts[[numHiddenLayers+1]]$holdA <- TRUE
    bestPar <- optim( randomInit, fn=fn, gr=gr, control=list(fnscale=-1, maxit=100, trace=TRUE, REPORT=1), method="BFGS",
                      u=uTotal, y=c(yInit,y), lu=lu, ly=ly, lag=lag, parSkeleton=optimInit, optimOpts=optimOpts, eps=eps, verbose=verbose )$par 
    
    optimOpts[[numHiddenLayers+1]]$holdA <- FALSE
    optimOpts <- lapply( optimOpts, "[[<-", "holdFirst", FALSE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdEps", FALSE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdSigma", FALSE )
    bestPar <- optim( bestPar, fn=fn, gr=gr, control=list(fnscale=-1, maxit=1000, trace=TRUE, REPORT=1), method="BFGS",
                      u=uTotal, y=c(yInit,y), lu=lu, ly=ly, lag=lag, parSkeleton=optimInit, optimOpts=optimOpts, eps=eps, verbose=verbose )$par
    
  } else{
    
    optimOpts <- lapply( optimOpts, "[[<-", "holdSigma", TRUE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdFirst", TRUE )
    bestPar <- optim( randomInit, fn=fn, gr=gr, control=list(fnscale=-1, maxit=100, trace=TRUE, REPORT=1), method="BFGS",
                      u=uTotal, y=c(yInit,y), lu=lu, ly=ly, lag=lag, parSkeleton=optimInit, optimOpts=optimOpts, eps=eps, verbose=verbose )$par 
    
    optimOpts <- lapply( optimOpts, "[[<-", "holdSigma", FALSE )
    optimOpts <- lapply( optimOpts, "[[<-", "holdFirst", FALSE )
    bestPar <- optim( bestPar, fn=fn, gr=gr, control=list(fnscale=-1, maxit=1000, trace=TRUE, REPORT=1), method="BFGS",
                      u=uTotal, y=c(yInit,y), lu=lu, ly=ly, lag=lag, parSkeleton=optimInit, optimOpts=optimOpts, eps=eps, verbose=verbose )$par
  }
  
  fn(bestPar,u=uTotal,y=c(yInit,y),lu,ly,lag,optimInit,optimOpts,eps=eps)
  
  cat("\n")
  plot( c(yInit, y), type="l", ylab="" )
  colorVec <- c( palette()[2:(numHiddenLayers+1)], "black" )
  for( h in 1:(numHiddenLayers+1) ){
    cat("Layer", h, "\n")
    cat("parF", model[[h]]$parF, "\n")
    cat("eps", model[[h]]$eps, "\n")
    if( h <= numHiddenLayers ){
      lines( model[[h]]$muVarPar, col=colorVec[h] )
      cat("sigma", model[[h]]$sigma, "\n")
      cat("muVarPar", range(model[[h]]$muVarPar), "\n")
      cat("lambdaVarPar", range(model[[h]]$lambdaVarPar), "\n")  
    }
    cat("\n")
  }
  legend( "topleft", legend=c( sapply( 1:numHiddenLayers, function(h) paste("Layer", h) ), "Training output" ), col=colorVec, lwd=2 )
  
  # Remove temp variables of the model
  for( h in 1:(numHiddenLayers+1) ){
    model[[h]]$auxZMat1 <- NULL
    model[[h]]$auxZMat2 <- NULL
    model[[h]]$auxZMat3 <- NULL
    model[[h]]$dPsi2ArrayCalc <- NULL
  }
  
  print( "Optimization done...", quote=FALSE )
  
  return( model )
}

gpRegression <- function( u=NULL, y=NULL, x=NULL, lu, ly, lag=1, M=0, numHiddenLayers=1, par=NULL, kernel="SE", robustFlag=FALSE, verbose=FALSE, ... ){
  
  kernel <- toupper( kernel )
  generateRandomParams <<- generateRandomParams_SE
  calcKernelMatrix <<- calcKernelMatrix_SE  
  calcGradLowerBound <<- calcGradLowerBound_SE
  numParDynKernel <- 1
  
  if( M <= 0 ){
    M <- round(0.1*length(y))
    print( paste("Using", M, "inducing points."), quote=FALSE )
  }  
  
  if( is.null(par) ){    
    if( verbose ){ print("Optimizing model...", quote=FALSE) }
    model <- optimizeHyperParameters( u=as.matrix(u), y=as.matrix(y), lu=lu, ly=ly, lag=lag, M=M, numHiddenLayers=numHiddenLayers, numParDynKernel=numParDynKernel, kernel=kernel, robustFlag=robustFlag, verbose=verbose )  
  }
  
  if( verbose ){ print("Estimation concluded.", quote=FALSE) }
  list( model=model, dataU=as.matrix(u), dataY=as.matrix(y), kernel=kernel, M=M, lu=lu, ly=ly, lag=lag, robustFlag=robustFlag )
}

# Not implemented
gpPredict <- cmpfun( function( m, u=NULL, y=NULL, x=NULL, verbose=TRUE, ... ){
  
  kernel <- toupper( m$kernel )
  numParDynKernel <- 1
  calcKernelMatrix <- calcKernelMatrix_SE  
  
  lu <- m$lu
  ly <- m$ly
  
  if( verbose ){ print("Formatting test data...", quote=FALSE) }
  if( is.null(x) ){
    formatted <- formatData( u, y, lu, ly )
    x <- formatted$x
    y <- formatted$y
  } 
  
  yEst <- rep(0, NROW(x))
  yVar <- rep(0, NROW(x))
  
  list( est=yEst, var=yVar )
} )

gpSimulation <- cmpfun( function( m, u=NULL, y=NULL, continueSimulation=FALSE, verbose=TRUE, ... ){
  
  kernel <- toupper( m$kernel )
  numParDynKernel <- 1
  calcKernelMatrix <- calcKernelMatrix_SE  
  
  if( verbose ){ print("Formatting test data...", quote=FALSE) }
  
  lu <- m$lu
  ly <- m$ly
  lag <- m$lag
  if( is.null(lag) ){
    lag <- 1
  }
  u <- as.matrix(u)
  
  numHiddenLayers <- length(m$model) - 1
  dataY <- m$dataY[-(1:ly)]
  N <- NROW(dataY)
  M <- m$M
  
  model <- m$model
  
  if( verbose ){ print("Performing free simulation...", quote=FALSE) }
  
  l <- max(lu,ly)
  simulation <- vector( "list", numHiddenLayers+1 )
  if( lu == 0 ){
    for( h in 1:numHiddenLayers ){
      simulation[[h]]$est <- rep(0, NROW(u)+ly)
      simulation[[h]]$var <- rep(model[[h]]$sigma, NROW(u)+ly)
      simulation[[h]]$est[1:ly] <- model[[h]]$muVarPar[(N+1):(N+ly)]  
      simulation[[h]]$var[1:ly] <- model[[h]]$lambdaVarPar[(N+1):(N+ly)]  
    }
    simulation[[numHiddenLayers+1]]$est <- c( dataY[(N+1):(N+l)], y[1:ly], rep(0, NROW(u)-ly) )
    simulation[[numHiddenLayers+1]]$var <- rep(model[[numHiddenLayers+1]]$sigma, NROW(u)+ly)
  } else{
    if( continueSimulation == FALSE ){
      for( h in 1:numHiddenLayers ){
        simulation[[h]]$est <- c( rep(0,l), rep(0, NROW(u)) )
        simulation[[h]]$var <- rep(model[[h]]$sigma, NROW(u)+l)
      }  
      simulation[[numHiddenLayers+1]]$est <- c( rep(0,l), y[1:ly], rep(0, NROW(u)-l) )
      simulation[[numHiddenLayers+1]]$var <- rep(model[[numHiddenLayers+1]]$sigma, NROW(u)+l)
      u <- rbind( as.matrix(rep(0,l)), u ) 
    } else{
      for( h in 1:numHiddenLayers ){
        simulation[[h]]$est <- rep(0, NROW(u)+l) 
        simulation[[h]]$var <- rep(model[[h]]$sigma, NROW(u)+l)
        simulation[[h]]$est[1:l] <- model[[h]]$muVarPar[(N+1):(N+l)]  
        simulation[[h]]$var[1:l] <- model[[h]]$lambdaVarPar[(N+1):(N+l)]  
      }
      simulation[[numHiddenLayers+1]]$est <- c( dataY[(N+1):(N+l)], y[1:l], rep(0, NROW(u)-l) )
      simulation[[numHiddenLayers+1]]$var <- rep(model[[numHiddenLayers+1]]$sigma, NROW(u)+l)
      u <- rbind( m$dataU[(N+1):(N+l),,drop=FALSE], u )  
    }
  }
  
  for( h in 1:numHiddenLayers ){
    model[[h]]$B <- model[[h]]$auxMatInv %*% t(model[[h]]$psi1) %*% model[[h]]$muVarPar[(ly+1):(N+ly)]
  }
  model[[numHiddenLayers+1]]$B <- model[[numHiddenLayers+1]]$auxMatInv %*% t(model[[numHiddenLayers+1]]$psi1) %*% dataY
  
  for( i in (l+1):length(simulation[[numHiddenLayers+1]]$est) ){
    for( h in 1:(numHiddenLayers+1) ){
      layer <- model[[h]]
      Q <- length(layer$parF[-1])
      iniX <- i-1; endX <- i-ly;
      
      if( h == 1 ){
        indu <- (ly-lu) + (i-lag+lu-l):(i-lag+1-l)
        if( lu > 0 ){
          muNewMat <- c( simulation[[h]]$est[iniX:endX], as.numeric(u[indu,]) )
        } else{
          muNewMat <- simulation[[h]]$est[iniX:endX]
        }
        sNewMat <- c( simulation[[h]]$var[iniX:endX], rep(0,lu*NCOL(u)) )
      } else if( h == numHiddenLayers+1 ){
        muNewMat <- simulation[[h-1]]$est[ 1+(iniX:endX) ]
        sNewMat <- simulation[[h-1]]$var[ 1+(iniX:endX) ] 
      } else{
        muNewMat <- c( simulation[[h]]$est[iniX:endX], simulation[[h-1]]$est[ 1+(iniX:endX) ])
        sNewMat <- c( simulation[[h]]$var[iniX:endX], simulation[[h-1]]$var[ 1+(iniX:endX) ] ) 
      }
      
      psiNew0 <- layer$sigmaF
      
      psiNew1 <- matrix( 0, nrow=1, ncol=M )
      n=1
      sMatDiag <- sNewMat
      psiNew1[n,] <- layer$sigmaF * apply( sapply( 1:Q, function(q)
        exp( -0.5*layer$w[q]*(muNewMat[q] - layer$zMat[q,])^2 / ( layer$w[q]*sMatDiag[q] + 1 ) ) / sqrt( layer$w[q]*sMatDiag[q] + 1 ) ), 1, prod ) 
      psiNew1 <- t(psiNew1)
      
      psiNew2 <- matrix( 0, nrow=M, ncol=M )
      n=1
      sMatDiag <- sNewMat
      for( m1 in 1:M ){
        psiNew2[m1,] <- (layer$sigmaF^2) * apply( sapply( 1:Q, function(q)
          exp( -0.25*layer$w[q]*( layer$zMat[q,m1] - layer$zMat[q,] )^2 - layer$w[q]*( muNewMat[q] - 0.5*(layer$zMat[q,m1] + layer$zMat[q,]) )^2 / ( 2*layer$w[q]*sMatDiag[q] + 1 ) ) / sqrt( 2*layer$w[q]*sMatDiag[q] + 1 ) ), 1, prod )
      }
      
      simulation[[h]]$est[i] <- as.numeric( t(layer$B) %*% psiNew1 )
      simulation[[h]]$var[i] <- ifelse( h < numHiddenLayers+1, model[[h]]$sigma, 0 ) + as.numeric(t(layer$B) %*% ( psiNew2 - psiNew1 %*% t(psiNew1) ) %*% layer$B) + psiNew0 - sum( (layer$KmmInv - layer$auxMatInv ) * t(psiNew2) ) 
    }
  }  
  
  if( m$robustFlag ){
    simulation[[numHiddenLayers+1]]$var <- simulation[[numHiddenLayers+1]]$var + median(model[[numHiddenLayers+1]]$bVarPar / model[[numHiddenLayers+1]]$aVarPar)
  } else{
    simulation[[numHiddenLayers+1]]$var <- simulation[[numHiddenLayers+1]]$var + simulation[[numHiddenLayers+1]]$sigma  
  }
  
  if( verbose ){ print("Free simulation concluded.", quote=FALSE) }
  
  return( list( est=simulation[[numHiddenLayers+1]]$est[-(1:l)], var=simulation[[numHiddenLayers+1]]$var[-(1:l)], simulation=simulation ) )  
} )

