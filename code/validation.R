library(ggplot2)
source("code/figuresFunctions.R", chdir=TRUE)

calculateMetrics <- function( true, est, var=NULL, verbose=TRUE ){
  true <- as.numeric(true)
  est <- as.numeric(est)
  err <- true - est 
  mse <- mean((err)^2)
  if( verbose ){ print(paste("MSE =", mse)) }
  rmse <- sqrt(mse)
  if( verbose ){ print(paste("RMSE =", rmse))}
  mrse <- sqrt( sum(err^2) / sum(true^2) )
  if( verbose ){ print(paste("MRSE =", mrse)) }
  nmse <- mean((err^2)/(mean(est)*mean(true)))
  if( verbose ){ print(paste("NMSE =", nmse)) }
  if( is.null(var) == FALSE ){
    smse <- mse / var( true )
    if( verbose ){ print(paste("SMSE =", smse)) }
    ld <- 0.5*log(2*pi) + (0.5/length(err)) * sum( log(var) + (err^2)/var )
    if( verbose ){ print(paste("NLD =", ld)) }
    list( err=err, mse=mse, rmse=rmse, mrse=mrse, smse=smse, ld=ld )
  } else{
    list( err=err, mse=mse, rmse=rmse, mrse=mrse, nmse=nmse )  
  }
}

plotValidation <- function( true, est, var=NULL, label="1 Step Prediction", figName="identification", sub=NULL, axisLabel=TRUE, font="timesnewroman", print=TRUE, save=TRUE ){
  # Fix warnings
  Step <- Output <- Signal <- Inf.Confidence.Interval <- Sup.Confidence.Interval <- NULL
  
  true <- as.numeric(true)
  est <- as.numeric(est)
  setFigureParams(font)
  if( is.null(var) ){
    var <- 0
  }
  t <- 1:length(true)
  df <- data.frame( Step=c(t,t),
                    Output=c(true, est),
                    Signal=c(rep("Real",length(est)),rep(label,length(est))),
                    Sup.Confidence.Interval=est+1.96*sqrt(var),
                    Inf.Confidence.Interval=est-1.96*sqrt(var))
  df$Signal <- factor(df$Signal, levels=unique(as.character(df$Signal)) )
  g <- ggplot( df, aes( x=Step, y=Output, col=Signal ) ) + geom_line(size=1) + geom_point(size=3) +
    scale_colour_manual(values=c("blue", "red"), 
                        labels=c("Real", label)) +
    geom_ribbon( aes(x=Step, ymin=Inf.Confidence.Interval, ymax=Sup.Confidence.Interval), df,
                 fill="red", col=NA, alpha=0.25 ) +
    guides(colour = guide_legend(override.aes = list(size = 7)))
  if( axisLabel == FALSE ){
    g <- g + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  }
  if( is.null(sub) == FALSE ){
    g <- g + ggtitle(sub)
  }
  if( save ){
    saveFigure( g, figName )  
  }
  if( print ){
    print( g )  
  } else{
    return( g )
  }
  
}
