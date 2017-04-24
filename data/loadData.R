library(R.matlab)

dataPath <- "data/"

loadData <- function( name ){
  
  out <- NULL
  
  if( name == "identificationExample" ){  
    identificationExample <- readMat( paste( dataPath, name, ".mat", sep="" ) )
    out <- list( dataName="Artificial example",
                 u=identificationExample$u, y=identificationExample$y, n=identificationExample$n, lu=1, ly=1 )
    
  } else if( length(grep( "identificationExample_30out", name )) > 0 ){  
    identificationExample <- readMat( paste( dataPath, name, ".mat", sep="" ) )
    out <- list( dataName="Artificial example 30% outliers",
                 u=identificationExample$u, y=identificationExample$y, n=identificationExample$n, lu=1, ly=1 )
    
  }                           
  
  out
}
