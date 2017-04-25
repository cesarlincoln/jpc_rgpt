library(Matrix)
library(compiler)

# Linear kernel

generateRandomParams_LIN <- function( numAtt ){  
  #   randomInit <- log( c( rep( 1/numAtt, numAtt ), rep( exp(-2), 2 ) ) )
  randomInit <- log( c( runif(numAtt+1, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_LIN <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_l <- par[1:D]
  sigma_c <- par[length(par)]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    .rowSums( t( t( x1[i,,drop=FALSE] * x2[j,,drop=FALSE] ) * sigma_l ), d1d2, D ) + sigma_c
  })
} )

calcRegularizationPenalty_LIN <- function( par, C ){
  C[1] * sum(par) + C[2] * sum(par*par) 
}

calcGradLogLikelihood_LIN <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  yLen <- length(y)
  sol <- exp(par)
  sigma_l <- sol[1:ncol(x)]
  sigma_c <- sol[ncol(x)+1]
  sigma_n <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  
  grSigmal <- rep(0, length(sigma_l))
  for( i in 1:length(sigma_l) ){
    grSigmal[i] <- 0.5 * sum( beta * t( x[,i,drop=FALSE] %*% t(x[,i,drop=FALSE]) ) )
    if( is.null( C ) == FALSE ){
      grSigmal[i] <- grSigmal[i] - C[1]*yLen - 2*C[2]*yLen*sol[1+i]  
    }
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmal, grSigmac, grSigman ) * sol
}

# Squared exponential kernel

generateRandomParams_SE <- function( numAtt ){  
  #   numAtt <- ncol(x)
  # randomInit <- log( c( 1, rep( 1/numAtt, numAtt ), exp(-2) ) )
  randomInit <- log( c( runif(1, min=0, max=1), runif(numAtt, min=0, max=1), runif(1, min=0, max=1) ) )
  
  #   R <- quantile( y, 0.95 )
  #   d <- sapply( 1:numAtt, function(i) median( as.matrix( dist( x[,i] ) ) ) )
  #   mu <- 0.5*log(R*d)
  #   sd <- 0.25*log(R/d)
  #   randomInit <- c( rnorm(1, mean=0, sd=1),
  #                    rnorm(numAtt, mean=mu, sd=sd),
  #                    rnorm(1, mean=-3, sd=1) )
}

calcKernelMatrix_SE <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  if( is.null(par) ){
    par <- c( 1, rep(1,D) )
  }
  sigma_f <- par[1]
  l <- par[2:length(par)]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    sigma_f * exp( -0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) )
  })
} )

calcRegularizationPenalty_SE <- calcRegularizationPenalty_LIN

calcGradLogLikelihood_SE <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  yLen <- length(y)
  sol <- exp(par)
  sigma <- sol[length(sol)] + eps
  sol <- sol[1:(length(sol)-1)]
  beta <- alpha %*% t(alpha) - KInv
  K2 <- K - sigma*diag(nrow(x))
  
  dSigmaf <- K2 / sol[1]
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sol[1]
  }
  
  grSigmaw <- rep(0, length(sol)-1)
  for( i in 1:(length(sol)-1) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*sol[1+i]
    }
  }
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigman ) * c(sol,sigma)
}

# Quadratic exponencial kernel with linear (single weight) and constant components

generateRandomParams_SELIN <- function( numAtt ){
  #   numAtt <- ncol(x)  
  #   randomInit <- log( c( 1,
  #                         rep( 1/numAtt, numAtt ),
  #                         rep( exp(-2), 2 ), 
  #                         exp(-2) ) )
  randomInit <- log( c( runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1),
                        runif(2, min=0, max=1), 
                        runif(1, min=0, max=1) ) )
  
  #   R <- quantile( y, 0.95 )
  #   d <- sapply( 1:numAtt, function(i) median( as.matrix( dist( x[,i] ) ) ) )
  #   mu <- 0.5*log(R*d)
  #   sd <- 0.25*log(R/d)
  #   randomInit <- c( rnorm(1, mean=0, sd=1),
  #                    rnorm(numAtt, mean=mu, sd=sd),
  #                    rnorm(2, mean=-2, sd=1), 
  #                    rnorm(1, mean=-3, sd=1) )
}

calcKernelMatrix_SELIN <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(D+1)]
  sigma_l <- par[D+2]
  sigma_c <- par[D+3]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    sigma_f * exp( -0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) ) + sigma_l * .rowSums( x1[i,,drop=FALSE] * x2[j,,drop=FALSE], d1d2, D ) + sigma_c
  })
} )

calcRegularizationPenalty_SELIN <- function( par, C ){
  C[1] * sum(par[1:(length(par)-1)]) + C[2] * sum(par[1:(length(par)-1)]*par[1:(length(par)-1)])  
}

calcGradLogLikelihood_SELIN <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){
  yLen <- length(y)
  sol <- exp(par)
  sigma_f <- sol[1]
  l <- sol[2:(ncol(x)+1)]
  sigma_l <- sol[ncol(x)+2]
  sigma_c <- sol[ncol(x)+3]
  sigma_n <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  xxMat <- x %*% t(x)
  K2 <- K - sigma_l * xxMat - sigma_c - sigma_n*diag(nrow(x))
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f
  }
  
  grSigmaw <- rep(0, length(l))  
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i]
    }
  }
  
  dSigmal <- xxMat
  grSigmal <- 0.5 * sum( beta * t(dSigmal) )
  if( is.null(C) == FALSE ){
    grSigmal <- grSigmal - C[1]*yLen - 2*C[2]*yLen*sigma_l
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmal, grSigmac, grSigman ) * sol
}

# Quadratic exponencial kernel with linear (multiple weights) and constant components 

generateRandomParams_SELINMULT <- function( numAtt ){  
  randomInit <- log( c( runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1), 
                        runif(numAtt, min=0, max=1), 
                        runif(1, min=0, max=1),
                        runif(1, min=0, max=1) ) )
}

calcKernelMatrix_SELINMULT <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(D+1)]
  sigma_l <- par[(D+2):(2*D+1)]
  sigma_c <- par[2*D+2]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    sigma_f * exp( -0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) ) + .rowSums( t( t( x1[i,,drop=FALSE] * x2[j,,drop=FALSE] ) * sigma_l ), d1d2, D ) + sigma_c
  })
} )

calcRegularizationPenalty_SELINMULT <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_SELINMULT <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){
  yLen <- length(y)
  sol <- exp(par)
  sigma_f <- sol[1]
  l <- sol[2:(ncol(x)+1)]
  sigma_l <- sol[(ncol(x)+2):(2*ncol(x)+1)]
  sigma_c <- sol[2*ncol(x)+2]
  sigma_n <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  xxMat <- ( x %*% diag(sigma_l) ) %*% t(x)    
  K2 <- K - xxMat - sigma_c - sigma_n*diag(nrow(x))
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f
  }
  
  grSigmaw <- rep(0, length(l))
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i]
    }
  }
  
  grSigmal <- rep(0, length(sigma_l))
  for( i in 1:length(sigma_l) ){
    grSigmal[i] <- 0.5 * sum( beta * t( x[,i,drop=FALSE] %*% t(x[,i,drop=FALSE]) ) )
    if( is.null(C) == FALSE ){
      grSigmal[i] <- grSigmal[i] - C[1]*yLen - 2*C[2]*yLen*sigma_l[i]
    }
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmal, grSigmac, grSigman ) * sol
}

# Quadratic exponencial kernel with periodic (unit period), linear (multiple weights) and constant components 

generateRandomParams_PER <- function( numAtt ){
  randomInit <- log( c( runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1), 
                        runif(numAtt, min=0, max=1), 
                        runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1),
                        runif(1, min=0, max=1),
                        runif(1, min=0, max=1) ) )
}

calcKernelMatrix_PER <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(D+1)]
  sigma_l <- par[(D+2):(2*D+1)]
  sigma_s <- par[2*D+2]
  ws <- par[(2*D+3):(3*D+2)]
  sigma_c <- par[3*D+3]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){
    sigma_f * exp( -0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) ) +
      sigma_s * exp( -2 * .rowSums( t( t( sin( 0.5*(x1[i,,drop=FALSE] - x2[j,,drop=FALSE]) )^2 ) * ws ), d1d2, D ) ) +
      .rowSums( t( t( x1[i,,drop=FALSE] * x2[j,,drop=FALSE] ) * sigma_l ), d1d2, D ) + sigma_c
  })
} )

calcRegularizationPenalty_PER <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_PER <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){
  D <- ncol(x)
  yLen <- length(y)
  sol <- exp(par)
  sigma_f <- sol[1]
  l <- sol[2:(D+1)]
  sigma_l <- sol[(D+2):(2*D+1)]
  sigma_s <- sol[2*D+2]
  ws <- sol[(2*D+3):(3*D+2)]
  sigma_c <- sol[3*D+3]
  sigma_n <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  xxMat <- ( x %*% diag(sigma_l) ) %*% t(x)
  d1d2 <- nrow(x)*nrow(x)  
  xExpSin <- outer( X=1:nrow(x), Y=1:nrow(x), FUN=function(i,j){
    sigma_s * exp( -2 * .rowSums( t( t( sin( 0.5*(x[i,,drop=FALSE] - x[j,,drop=FALSE]) )^2 ) * ws ), d1d2, D ) )
  })
  
  K2 <- K - xxMat - xExpSin - sigma_c - sigma_n*diag(nrow(x))
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f  
  }
  
  grSigmaw <- rep(0, length(l))
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i]  
    }
  }
  
  grSigmal <- rep(0, length(sigma_l))
  for( i in 1:length(sigma_l) ){
    grSigmal[i] <- 0.5 * sum( beta * t( x[,i,drop=FALSE] %*% t(x[,i,drop=FALSE]) ) )
    if( is.null(C) == FALSE ){
      grSigmal[i] <- grSigmal[i] - C[1]*yLen - 2*C[2]*yLen*sigma_l[i]  
    }
  }
  
  dSigmas <- xExpSin / sigma_s
  grSigmas <- 0.5 * sum( beta * t(dSigmas) )
  if( is.null(C) == FALSE ){
    grSigmas <- grSigmas - C[1]*yLen - 2*C[2]*yLen*sigma_s
  }
  
  grSigmaWs <- rep(0, length(ws))
  for( i in 1:length(ws) ){
    grSigmaWs[i] <- 0.5 * sum( beta * t( -2 * ( sin( 0.5 * as.matrix( outer( x[,i], x[,i], function(p,q) abs(p-q) ) ) )^2 ) * xExpSin ) )
    if( is.null(C) == FALSE ){
      grSigmaWs[i] <- grSigmaWs[i] - C[1]*yLen - 2*C[2]*yLen*ws[i]  
    }
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmal, grSigmas, grSigmaWs, grSigmac, grSigman ) * sol
}

# Quadratic exponencial kernel with periodic (variable period), linear (multiple weights) and constant components

generateRandomParams_PERVAR <- function( numAtt ){
  randomInit <- log( c( runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1), 
                        runif(numAtt, min=0, max=1), 
                        runif(1, min=0, max=1),
                        runif(numAtt, min=0, max=1),
                        runif(1, min=0, max=2),
                        runif(1, min=0, max=1),
                        runif(1, min=0, max=1) ) )
}

calcKernelMatrix_PERVAR <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(D+1)]
  sigma_l <- par[(D+2):(2*D+1)]
  sigma_s <- par[2*D+2]
  ws <- par[(2*D+3):(3*D+2)]
  ps <- par[3*D+3]
  sigma_c <- par[3*D+4]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){
    sigma_f * exp( -0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) ) +
      sigma_s * exp( -2 * .rowSums( t( t( sin( pi*ps*(x1[i,,drop=FALSE] - x2[j,,drop=FALSE]) )^2 ) * ws ), d1d2, D ) ) +
      .rowSums( t( t( x1[i,,drop=FALSE] * x2[j,,drop=FALSE] ) * sigma_l ), d1d2, D ) + sigma_c
  })
} )

calcRegularizationPenalty_PERVAR <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_PERVAR <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){
  D <- ncol(x)
  yLen <- length(y)
  sol <- exp(par)
  sigma_f <- sol[1]
  l <- sol[2:(D+1)]
  sigma_l <- sol[(D+2):(2*D+1)]
  sigma_s <- sol[2*D+2]
  ws <- sol[(2*D+3):(3*D+2)]
  ps <- sol[3*D+3]
  sigma_c <- sol[3*D+4]
  sigma_n <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  xxMat <- ( x %*% diag(sigma_l) ) %*% t(x)
  
  d1d2 <- nrow(x)*nrow(x)  
  xExpSin <- outer( X=1:nrow(x), Y=1:nrow(x), FUN=function(i,j){
    sigma_s * exp( -2 * .rowSums( t( t( sin( pi*ps*(x[i,,drop=FALSE] - x[j,,drop=FALSE]) )^2 ) * ws ), d1d2, D ) )
  })
  dPs <- xExpSin * outer( X=1:nrow(x), Y=1:nrow(x), FUN=function(i,j){
    -4 * .rowSums( t( t( 0.5 * sin( 2*pi*ps*(x[i,,drop=FALSE] - x[j,,drop=FALSE]) ) * pi * (x[i,,drop=FALSE] - x[j,,drop=FALSE]) ) * ws ), d1d2, D )
  })
  
  K2 <- K - xxMat - xExpSin - sigma_c - sigma_n*diag(nrow(x))
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f 
  }
  
  grSigmaw <- rep(0, length(l))
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i] 
    }
  }
  
  grSigmal <- rep(0, length(sigma_l))
  for( i in 1:length(sigma_l) ){
    grSigmal[i] <- 0.5 * sum( beta * t( x[,i,drop=FALSE] %*% t(x[,i,drop=FALSE]) ) )
    if( is.null(C) == FALSE ){
      grSigmal[i] <- grSigmal[i] - C[1]*yLen - 2*C[2]*yLen*sigma_l[i] 
    }
  }
  
  dSigmas <- xExpSin / sigma_s
  grSigmas <- 0.5 * sum( beta * t(dSigmas) )
  if( is.null(C) == FALSE ){
    grSigmas <- grSigmas - C[1]*yLen - 2*C[2]*yLen*sigma_s
  }
  
  grSigmaWs <- rep(0, length(ws))
  for( i in 1:length(ws) ){
    grSigmaWs[i] <- 0.5 * sum( beta * t( -2 * ( sin( pi * ps * as.matrix( outer( x[,i], x[,i], function(p,q) abs(p-q) ) ) )^2 ) * xExpSin ) )
    if( is.null(C) == FALSE ){
      grSigmaWs[i] <- grSigmaWs[i] - C[1]*yLen - 2*C[2]*yLen*ws[i] 
    }
  }
  
  grPs <- 0.5 * sum( beta * t(dPs) )
  if( is.null(C) == FALSE ){
    grPs <- grPs - C[1]*yLen - 2*C[2]*yLen*ps
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmal, grSigmas, grSigmaWs, grPs, grSigmac, grSigman ) * sol
}

# Matern Kernel with nu=3/2

generateRandomParams_MATERN <- function( numAtt ){
  randomInit <- log( c( runif(1, min=0, max=1), runif(numAtt+1, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_MATERN <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(1+D)]
  sigma_c <- par[length(par)]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    aux <- sqrt( 3 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) )
    sigma_f * exp( - aux ) * ( 1 + aux ) + sigma_c
  })
} )

calcRegularizationPenalty_MATERN <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_MATERN <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  yLen <- length(y)
  sol <- exp(par)    
  sigma <- sol[length(sol)]
  sigma_f <- sol[1]
  l <- sol[2:(length(sol)-2)]
  sigma_c <- sol[length(sol)-1]
  
  beta <- alpha %*% t(alpha) - KInv
  K2 <- K - sigma_c - sigma*diag(nrow(x))
  
  D <- ncol(x)
  d1d2 <- nrow(x)*nrow(x)  
  K3 <- outer( X=1:nrow(x), Y=1:nrow(x), FUN=function(i,j){ 
    sqrt(3) * sigma_f * exp( - sqrt( 3 * .rowSums( t( t( (x[i,,drop=FALSE] - x[j,,drop=FALSE])^2 ) * l ), d1d2, D ) ) )
  })
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f
  }
  
  grSigmaw <- rep(0, length(l))
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * sqrt(3) * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K3 ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i]
    }
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmac, grSigman ) * sol
}

# Rational Quadratic Kernel

generateRandomParams_RQ <- function( numAtt ){  
  randomInit <- log( c( runif(1, min=0, max=1), runif(numAtt+1, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_RQ <- cmpfun( function( x1, x2, par ){  
  D <- ncol(x1)
  sigma_f <- par[1]
  l <- par[2:(D+1)]
  a <- par[length(par)]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    sigma_f * ( 1 + 0.5 * .rowSums( t( t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 ) * l ), d1d2, D ) / a )^(-a)
  })
  
  #   d1 <- nrow(x1)
  #   d2 <- nrow(x2)
  #   res <- matrix(-1, nrow=d1, ncol=d2)  
  #   if( d1 == d2 ){
  #     for( i in seq_len(d1) ){
  #       xi <- x1[i,]
  #       for( j in i:d2 ){
  #         res[i,j] <- sigma_f * ( 1 + 0.5*sum( l * (xi - x2[j,])^2) / a )^(-a) 
  #       }
  #     }
  #     res <- as.matrix(forceSymmetric(res))
  #   } else{
  #     for( i in seq_len(d1) ){
  #       xi <- x1[i,] 
  #       for( j in seq_len(d2) ){
  #         res[i,j] <- sigma_f * ( 1 + 0.5*sum( l * (xi - x2[j,])^2) / a )^(-a) 
  #       }
  #     }  
  #   }  
  #   
  #   res
} )

calcRegularizationPenalty_RQ <- function( par, C ){
  aux <- par[1:(length(par)-1)]
  C[1] * sum(aux) + C[2] * sum(aux*aux) 
}

calcGradLogLikelihood_RQ <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  D <- ncol(x)
  yLen <- length(y)
  sol <- exp(par)
  sigma_f <- sol[1]
  l <- sol[2:(D+1)]
  a <- sol[length(sol)-1]
  sigman <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  K2 <- K - sigman*diag(nrow(x))
  
  d1d2 <- nrow(x)*nrow(x)  
  Z <- outer( X=1:nrow(x), Y=1:nrow(x), FUN=function(i,j){ 
    1 + 0.5 * .rowSums( t( t( (x[i,,drop=FALSE] - x[j,,drop=FALSE])^2 ) * l ), d1d2, D ) / a
  })
  
  dSigmaf <- K2 / sigma_f
  grSigmaf <- 0.5 * sum( beta * t(dSigmaf) )
  if( is.null(C) == FALSE ){
    grSigmaf <- grSigmaf - C[1]*yLen - 2*C[2]*yLen*sigma_f
  }
  
  grSigmaw <- rep(0, length(l))
  for( i in 1:length(l) ){
    grSigmaw[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2 / Z ) )
    if( is.null(C) == FALSE ){
      grSigmaw[i] <- grSigmaw[i] - C[1]*yLen - 2*C[2]*yLen*l[i]
    }
  }
  
  dSigmaa <- K2 * ( 1 - 1/Z - log(Z) )
  grSigmaa <- 0.5 * sum( beta * t(dSigmaa) )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf, grSigmaw, grSigmaa, grSigman ) * sol
}

# Spectral Mixture Kernel

generateRandomParams_SM <- function( numAtt ){
  Q <- 5
  randomInit <- log( c( runif(Q, min=0, max=1), runif(2*Q*numAtt, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_SM <- cmpfun( function( x1, x2, par ){  
  Q <- 5
  
  numAtt <- ncol(x1)
  w <- par[1:Q]
  tmu2pi <- 2*pi*t( matrix( par[(Q+1):(Q+Q*numAtt)], nrow=Q, ncol=numAtt ) )
  tv2pi2Negative <- -2*pi*pi*t( matrix( par[(Q+Q*numAtt+1):length(par)], nrow=Q, ncol=numAtt ) )
  
  onesVec <- rep( 1, Q )
  d1 <- nrow(x1)
  d2 <- nrow(x2)
  res <- matrix(-1, nrow=d1, ncol=d2)  
  if( d1 == d2 ){
    for( i in seq_len(d1) ){
      xi <- x1[i,]
      for( j in i:d2 ){
        r <- xi - x2[j,]
        prodAux <- exp(t(r*r*tv2pi2Negative))
        colAux <- onesVec
        for( k in 1:numAtt ){ colAux <- prodAux[,k] * colAux }
        res[i,j] <- sum( w * cos(t(r*tmu2pi)) * colAux )
      }
    }
    res <- as.matrix(forceSymmetric(res))
  } else{
    for( i in seq_len(d1) ){
      xi <- x1[i,] 
      for( j in seq_len(d2) ){
        r <- xi - x2[j,]
        prodAux <- exp(t(r*r*tv2pi2Negative))
        colAux <- onesVec
        for( k in 1:numAtt ){ colAux <- prodAux[,k] * colAux }
        res[i,j] <- sum( w * cos(t(r*tmu2pi)) * colAux )
      }
    }  
  }  
  
  res
} )

calcRegularizationPenalty_SM <- function( par, C ){
  0
}

# This gradient is wrong because of a typo in the original paper...
calcGradLogLikelihood_SM <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  numAtt <- ncol(x)
  yLen <- length(y)
  sol <- exp(par)
  Q <- 5
  
  w <- sol[1:Q]
  tmu <- t( matrix( sol[(Q+1):(Q+Q*numAtt)], nrow=Q, ncol=numAtt ) )
  tv <- t( matrix( sol[(Q+Q*numAtt+1):(length(sol)-1)], nrow=Q, ncol=numAtt ) )
  sigman <- sol[length(sol)] + eps
  
  beta <- alpha %*% t(alpha) - KInv
  
  zeroMat <- matrix( 0, nrow=nrow(K), ncol=ncol(K) )  
  dwMatList <- vector( "list", Q )
  dwMatList <- lapply( dwMatList, function(i) zeroMat )
  dmuList <- vector( "list", Q*ncol(x) )  
  dmuList <- lapply( dmuList, function(i) zeroMat )
  dvList <- dmuList
  
  tmu2pi <- 2*pi*tmu
  tv2pi2Negative <- -2*pi*pi*tv
  
  onesVec <- rep( 1, Q )
  indexMat <- matrix( 1:(Q*numAtt), nrow=Q, ncol=numAtt )
  for( i in seq_len(nrow(x)) ){    
    xi <- x[i,]
    for( j in i:nrow(x) ){
      r <- xi - x[j,]
      expAux <- exp(t(r*r*tv2pi2Negative))
      aux <- t(r*tmu2pi)
      cosAux <- cos(aux)
      sinAux <- sin(aux)
      
      colAux  <- onesVec
      prodAux <- expAux
      for( k in 1:numAtt ){ colAux <- prodAux[,k] * colAux }
      prodAux <- colAux
      
      for( q in 1:Q ){
        dwMatList[[q]][i,j] <- prodAux[q] * cosAux
      }
      for( q in 1:nrow(indexMat) ){
        prodAuxQ <- prodAux[q]
        for( p in 1:ncol(indexMat) ){
          aux <- 2*pi * w[q] * cosAux[q] * r[p]
          dmuList[[indexMat[q,p]]][i,j] <- - aux * ( prodAuxQ / cosAux[q,p] ) * sinAux[q,p]
          dvList[[indexMat[q,p]]][i,j] <- - pi * r[p] * aux * prodAuxQ
        }
      }
    }
  }
  dwMatList <- lapply( dwMatList, function(k) as.matrix(forceSymmetric(k)) )
  dmuList <- lapply( dmuList, function(k) as.matrix(forceSymmetric(k)) )
  dvList <- lapply( dvList, function(k) as.matrix(forceSymmetric(k)) )
  
  grW <- sapply( dwMatList, function(k) 0.5 * sum( beta * t(k) ) )  
  grMu <- sapply( dmuList, function(k) 0.5 * sum( beta * t(k) ) )  
  grV <- sapply( dvList, function(k) 0.5 * sum( beta * t(k) ) )  
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grW, grMu, grV, grSigman ) * sol
}

# Double squared exponential kernel

generateRandomParams_DSE <- function( numAtt ){  
  randomInit <- log( c( runif(1, min=0, max=1), runif(numAtt, min=0, max=1), runif(1, min=0, max=1), runif(numAtt, min=0, max=1), runif(2, min=0, max=1) ) )
}

calcKernelMatrix_DSE <- cmpfun( function( x1, x2, par=NULL ){  
  D <- ncol(x1)
  sigma_f1 <- par[1]
  l1 <- par[2:(1+D)]
  
  sigma_f2 <- par[2+D]
  l2 <- par[(3+D):(length(par)-1)]
  
  sigma_c <- par[length(par)]
  
  d1d2 <- nrow(x1)*nrow(x2)  
  outer( X=1:nrow(x1), Y=1:nrow(x2), FUN=function(i,j){ 
    aux <- t( (x1[i,,drop=FALSE] - x2[j,,drop=FALSE])^2 )
    sigma_f1 * exp( -0.5 * .rowSums( t( aux * l1 ), d1d2, D ) ) +
      sigma_f2 * exp( -0.5 * .rowSums( t( aux * l2 ), d1d2, D ) ) + sigma_c
  })
  
  #   d1 <- nrow(x1)
  #   d2 <- nrow(x2)
  #   res <- matrix(-1, nrow=d1, ncol=d2)  
  #   if( d1 == d2 ){
  #     for( i in seq_len(d1) ){
  #       xi <- x1[i,]
  #       for( j in i:d2 ){
  #         res[i,j] <- sigma_f1 * exp( -0.5 * sum( l1 * (xi - x2[j,])^2) ) + sigma_f2 * exp( -0.5 * sum( l2 * (xi - x2[j,])^2) ) + sigma_c
  #       }
  #     }
  #     res <- as.matrix(forceSymmetric(res))
  #   } else{
  #     for( i in seq_len(d1) ){
  #       xi <- x1[i,] 
  #       for( j in seq_len(d2) ){
  #         res[i,j] <- sigma_f1 * exp( -0.5 * sum( l1 * (xi - x2[j,])^2) ) + sigma_f2 * exp( -0.5 * sum( l2 * (xi - x2[j,])^2) ) + sigma_c
  #       }
  #     }  
  #   }  
  #   
  #   res
} )

calcRegularizationPenalty_DSE <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_DSE <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  yLen <- length(y)
  sol <- exp(par)
  sigma <- sol[length(sol)] + eps
  sol <- sol[1:(length(sol)-1)]
  
  sigma_f1 <- sol[1]
  l1 <- sol[2:(1+ncol(x))]  
  sigma_f2 <- sol[2+ncol(x)]
  l2 <- sol[(3+ncol(x)):(length(sol)-1)]  
  sigma_c <- sol[length(sol)]
  
  beta <- alpha %*% t(alpha) - KInv
  
  K2_1 <- calcKernelMatrix_SE( x, x, c(sigma_f1, l1) )
  K2_2 <- calcKernelMatrix_SE( x, x, c(sigma_f2, l2) )
  
  dSigmaf1 <- K2_1 / sigma_f1
  grSigmaf1 <- 0.5 * sum( beta * t(dSigmaf1) )
  if( is.null(C) == FALSE ){
    grSigmaf1 <- grSigmaf1 - C[1]*yLen - 2*C[2]*yLen*sigma_f1
  }
  
  grSigmaw1 <- rep(0, ncol(x))
  for( i in 1:length(grSigmaw1) ){
    grSigmaw1[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2_1 ) )
    if( is.null(C) == FALSE ){
      grSigmaw1[i] <- grSigmaw1[i] - C[1]*yLen - 2*C[2]*yLen*l1[i]
    }
  }
  
  dSigmaf2 <- K2_2 / sigma_f2
  grSigmaf2 <- 0.5 * sum( beta * t(dSigmaf2) )
  if( is.null(C) == FALSE ){
    grSigmaf2 <- grSigmaf2 - C[1]*yLen - 2*C[2]*yLen*sigma_f2
  }
  
  grSigmaw2 <- rep(0, ncol(x))
  for( i in 1:length(grSigmaw2) ){
    grSigmaw2[i] <- 0.5 * sum( beta * t( - 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * K2_2 ) )
    if( is.null(C) == FALSE ){
      grSigmaw2[i] <- grSigmaw2[i] - C[1]*yLen - 2*C[2]*yLen*l2[i]
    }
  }
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmaf1, grSigmaw1, grSigmaf2, grSigmaw2, grSigmac, grSigman ) * c(sol,sigma)
}

# Additive squared expontential kernel

generateRandomParams_ADD <- function( numAtt ){  
  R <- min(numAtt,10)
  randomInit <- log( c( runif(R, min=0, max=1), runif(numAtt, min=0, max=1), runif(numAtt, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_ADD <- cmpfun( function( x1, x2, par=NULL ){  
  numAtt <- ncol(x1)
  D <- numAtt
  R <- min(numAtt,10)
  
  sigma_r <- par[1:R]
  
  sigma_f <- par[(R+1):(R+numAtt)]
  l <- par[(R+numAtt+1):length(par)]
  
  d1 <- nrow(x1)
  d2 <- nrow(x2)
  d1d2 <- nrow(x1)*nrow(x2)  
  z <- vector( "list", numAtt )
  for( d in 1:numAtt ){    
    x1d <- x1[,d]; x2d <- x2[,d]
    sigma_fd <- sigma_f[d]; ld <- l[d]
    
    z[[d]] <- outer( X=1:length(x1d), Y=1:length(x2d), FUN=function(i,j){       
      sigma_fd * exp( -0.5 * (x1d[i] - x2d[j])^2 * ld )
    })
  }
  assign( "zList", z )
  
  e <- vector( "list", R+1 ); e[[1]] <- matrix( 1, nrow=nrow(x1), ncol=nrow(x2) )        
  for( n in 1:R ){   
    e[[n+1]] <- 0     
    for( k in 1:n ){
      e[[n+1]] <- e[[n+1]] + (-1)^(k-1) * e[[n-k+1]] * Reduce( "+", lapply(z, function(z) z^k ) )    
    }
    e[[n+1]] <- e[[n+1]] / n          
  }     
  assign( "eList", e )
  
  res <- 0   
  for( k in 1:R ){
    res <- res + sigma_r[k] * e[[k+1]]    
  }
  
  res
} )

calcRegularizationPenalty_ADD <- calcRegularizationPenalty_LIN

calcGradLogLikelihood_ADD <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  numAtt <- ncol(x)
  R <- min(numAtt,10)
  
  yLen <- length(y)
  sol <- exp(par)
  sigma <- sol[length(sol)] + eps
  sigma_r <- sol[1:R]  
  sigma_f <- sol[(R+1):(R+numAtt)]
  l <- sol[(R+numAtt+1):(length(sol)-1)]
  
  beta <- alpha %*% t(alpha) - KInv
  
  grSigmar <- sapply( 1:R, function(i) 0.5 * sum( beta * t(eList[[i+1]]) ) )
  
  grSigmaf <- rep(0, numAtt)
  grSigmaw <- rep(0, numAtt)
  oneMatrix <- matrix( 1, nrow=nrow(x), ncol=nrow(x) )
  for( i in 1:numAtt ){
    
    e <- vector( "list", R+1 ); e[[1]] <- oneMatrix; 
    for( n in 1:R ){   
      e[[n+1]] <- 0   
      
      for( k in 1:n ){ 
        e[[n+1]] <- e[[n+1]] + (-1)^(k-1) * e[[n-k+1]] * Reduce( "+", lapply( zList[-i], function(z) z^k ) )  
      }
      e[[n+1]] <- e[[n+1]] / n          
    }     
    
    auxGrad <- 0
    for( k in 1:R ){
      auxGrad <- auxGrad + sigma_r[k] * e[[k]]    
    }
    
    dSigmaf <- auxGrad * ( zList[[i]] / sigma_f[i] )
    grSigmaf[i] <- 0.5 * sum( beta * t(dSigmaf) )
    
    grSigmaw[i] <- 0.5 * sum( beta * t( auxGrad * (- 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * zList[[i]] ) ) )
  }
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmar, grSigmaf, grSigmaw, grSigman ) * sol
}

# Additive squared expontential kernel with linear (single weight) and constant components

generateRandomParams_ADDLIN <- function( numAtt ){  
  R <- min(numAtt,10)
  randomInit <- log( c( runif(R, min=0, max=1), runif(numAtt, min=0, max=1), runif(numAtt, min=0, max=1), runif(2, min=0, max=1), runif(1, min=0, max=1) ) )
}

calcKernelMatrix_ADDLIN <- cmpfun( function( x1, x2, par=NULL ){  
  numAtt <- ncol(x1)
  D <- numAtt
  R <- min(numAtt,10)
  
  sigma_r <- par[1:R]
  
  sigma_f <- par[(R+1):(R+numAtt)]
  l <- par[(R+numAtt+1):(length(par)-2)]
  sigma_l <- par[length(par)-1]
  sigma_c <- par[length(par)]
  
  d1 <- nrow(x1)
  d2 <- nrow(x2)
  d1d2 <- nrow(x1)*nrow(x2)   
  z <- vector( "list", numAtt )
  for( d in 1:numAtt ){    
    x1d <- x1[,d]; x2d <- x2[,d]
    sigma_fd <- sigma_f[d]; ld <- l[d]
    
    z[[d]] <- outer( X=1:length(x1d), Y=1:length(x2d), FUN=function(i,j){       
      sigma_fd * exp( -0.5 * (x1d[i] - x2d[j])^2 * ld )
    })
  }
  
  assign( "zList", z )
  
  e <- vector( "list", R+1 ); e[[1]] <- matrix( 1, nrow=nrow(x1), ncol=nrow(x2) )        
  for( n in 1:R ){   
    e[[n+1]] <- 0     
    for( k in 1:n ){
      e[[n+1]] <- e[[n+1]] + (-1)^(k-1) * e[[n-k+1]] * Reduce( "+", lapply(z, function(z) z^k ) )    
    }
    e[[n+1]] <- e[[n+1]] / n          
  }     
  assign( "eList", e )
  
  res <- 0   
  for( k in 1:R ){
    res <- res + sigma_r[k] * e[[k+1]]    
  }
  
  res <- res + sigma_l * (x1 %*% t(x2)) + sigma_c
  
  res
} )

calcRegularizationPenalty_ADDLIN <- calcRegularizationPenalty_SELIN

calcGradLogLikelihood_ADDLIN <- function( par, x, y, eps=0, C=NULL, verbose=FALSE ){    
  numAtt <- ncol(x)
  R <- min(numAtt,10)
  
  yLen <- length(y)
  sol <- exp(par)
  sigma <- sol[length(sol)] + eps
  sigma_r <- sol[1:R]  
  sigma_f <- sol[(R+1):(R+numAtt)]
  l <- sol[(R+numAtt+1):(length(sol)-3)]
  sigma_l <- sol[length(sol)-2]
  sigma_c <- sol[length(sol)-1]
  
  beta <- alpha %*% t(alpha) - KInv
  
  grSigmar <- sapply( 1:R, function(i) 0.5 * sum( beta * t(eList[[i+1]]) ) )
  
  grSigmaf <- rep(0, numAtt)
  grSigmaw <- rep(0, numAtt)
  oneMatrix <- matrix( 1, nrow=nrow(x), ncol=nrow(x) )
  for( i in 1:numAtt ){
    
    e <- vector( "list", R+1 ); e[[1]] <- oneMatrix; 
    for( n in 1:R ){   
      e[[n+1]] <- 0   
      
      for( k in 1:n ){ 
        e[[n+1]] <- e[[n+1]] + (-1)^(k-1) * e[[n-k+1]] * Reduce( "+", lapply( zList[-i], function(z) z^k ) )  
      }
      e[[n+1]] <- e[[n+1]] / n          
    }     
    
    auxGrad <- 0
    for( k in 1:R ){
      auxGrad <- auxGrad + sigma_r[k] * e[[k]]    
    }
    
    dSigmaf <- auxGrad * ( zList[[i]] / sigma_f[i] )
    grSigmaf[i] <- 0.5 * sum( beta * t(dSigmaf) )
    
    grSigmaw[i] <- 0.5 * sum( beta * t( auxGrad * (- 0.5 * outer( x[,i], x[,i], function(p,q) (p-q)^2 ) * zList[[i]] ) ) )
  }
  
  dSigmal <- x %*% t(x)
  grSigmal <- 0.5 * sum( beta * t(dSigmal) )
  
  grSigmac <- 0.5 * sum( beta )
  
  dSigman <- diag(nrow(x))
  grSigman <- 0.5 * sum( beta * t(dSigman) )
  
  c( grSigmar, grSigmaf, grSigmaw, grSigmal, grSigmac, grSigman ) * sol
}

