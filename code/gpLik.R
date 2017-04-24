
calcLogLikelihood <- function( par, x, y, eps=0, verbose=FALSE ){    
  sol <- exp(par)
  #   if( verbose ){
  #     print(sol)
  #   }
  sigma <- sol[length(sol)] + eps
  sol <- sol[1:(length(sol)-1)]
  K <- calcKernelMatrix( x, x, sol ) + sigma*diag(nrow(x))
  cholK <- NA
  tryCatch( cholK <- chol(K), error=function(e) NA )
  # cholK <- jitchol(K)
  if( length(cholK) == 1 && is.na(cholK) ){ return(NA) }
  KInv <- chol2inv(cholK)
  
  K <<- K
  KInv <<- KInv
  alpha <- KInv %*% y
  alpha <<- alpha
  
  res <- -0.5 * t(y) %*% alpha - sum(log(diag(cholK))) - 0.5*length(y)*log(2*pi)  
  #   if( verbose ){
  #     print(res)
  #   }
  res
}