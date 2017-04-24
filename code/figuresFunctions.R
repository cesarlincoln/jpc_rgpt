library(grid)

setFigureParams <- function( font="timesnewroman", size=35 ){
  if(is.null(font)){
    theme_set( theme_bw( base_size = size ) ) # 25 35  
  } else{
    theme_set( theme_bw( base_size = size, base_family=font ) )  # 25 35  
  }
  theme_update( axis.title.y=element_text( angle=90, vjust=0.3 ),
                axis.title.x=element_text( vjust=0 ),
                plot.title = element_text( size=rel(0.8) ),
                legend.title=element_blank(), legend.key=element_rect(colour="white"), legend.position="top",
                legend.key.width=unit(2,"line"),
                legend.spacing.x = unit(0, "line"), legend.spacing.y = unit(0, "line"), plot.margin = unit(c(0,0.5,0.5,0.5), "lines") )
                # legend.margin = unit(0, "line"), plot.margin = unit(c(0,0.5,0.5,0.5), "lines") )
}

saveFigure <- function( fig, figName, plot=FALSE, main="", outDir=getwd(), width=1024, height=768, eps=TRUE ){  
  while( dev.cur() != 1 ){ dev.off() }
  figName <- paste( outDir, figName, sep="")
  png( filename=paste(figName,".png",sep=""), width, height )
  if( plot ){
    par( family="timesnewroman", cex.axis=2, cex.lab=2, cex.main=2 )
    plot(fig, main=main, lwd=4)  
  } else{
    print(fig)  
  }
  while( dev.cur() != 1 ){ dev.off() }
  if( eps ){
    system( paste( "sam2p ", figName, ".png EPS: ", figName, ".eps", sep="" ), ignore.stdout=TRUE, ignore.stderr=TRUE )
  }
}