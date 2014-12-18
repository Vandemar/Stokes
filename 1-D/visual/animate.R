#!/usr/bin/R

library(animation)

ts_ctl <- length(unique(data$ts))

pltData <- function() {
#for(i in 0:(ts_ctl-1)) {
for(i in 0:(100)) {
#  pdf(paste("Upwind_h256",i,".pdf",sep=""),height=4,width=6.5)
  points <- subset(data, ts==(ts_ctl-1))
  x <- points$x
  y <- points$u
  plot(x, y, type='l', ylim=c(0,1), main=) 
#  dev.off()
}
}

makePlot <- function(outputFile, inputFile ) {
  data <- read.table(inputFile, col.names=c("ts", "x", "u"), sep=',')
  oopt = ani.options(interval = 0, nmax = (ts_ctl-1))
  saveGIF(pltData(), movie.name = outputFile, interval = 0.1, width = 580, height = 400)
#  saveGIF(pltData(), interval = 0.1, width = 580, height = 400)
  ani.options(oopt)
}
