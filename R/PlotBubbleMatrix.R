PlotBubbleMatrix <-
function(Q, main="", special=Inf, cex=1){
diag(Q) <- 0
	Q<-Q/max(Q)

  plot(x=range(.5,.5+dim(Q)[2]),y=-range(.5, .5+dim(Q)[1]), xlab="", ylab="", type="n", main=main,xaxt='n',yaxt='n', asp=1,bty="n")
  	axis(side=2, at=-sequence(dim(Q)[1]), labels=rownames(Q), las=2, cex.axis=cex)
	axis(side=3, at=sequence(dim(Q)[2]), labels=colnames(Q), las=2, cex.axis=cex)

  for (i in sequence(dim(Q)[2])) {
  	for (j in sequence(dim(Q)[1])) {
  		bg="black"
  		if(i%in% special || j %in% special) {
  			bg="red"
  		}
  		if(Q[j,i]>0) {
  			symbols(x=i, y=-j, circles=sqrt(Q[j,i])/(2.1*sqrt(max(Q))), inches=FALSE, add=TRUE, fg=bg, bg=bg)
  		}
  	}	
  }
}
