library(expm)
library(miscTools)
library(igraph)

ComputePathwayProbability <- function(Q, from, to, correct.diag=TRUE, total.time=1000000) {
	
	if(correct.diag) {
		diag(Q) <- 0
		diag(Q) <- -rowSums(Q)
	}
	if(is.null(rownames(Q))) {
		rownames(Q) <- as.character(sequence(dim(Q)[1]))
		colnames(Q) <- rownames(Q)
	} 
	Q.expanded <- Q
	Q.expanded[ , to] <- 0
	Q.expanded[to , ] <- 0

	num.sources <- dim(Q)[1] - 1
	for (i in sequence(num.sources-1)) {
		Q.expanded <- insertRow(Q.expanded, to, 0)
		Q.expanded <- insertCol(Q.expanded, to, 0)
	}
	Q.expanded[ , to] <- 0
	Q.expanded[to , ] <- 0
	#ok, now we have a matrix with empty rows and cols for the absorbing states
	from.new <- from
	if (from > to) {
		from.new <- from + num.sources - 1 #to deal with expansion	
	}
	sources.original <- sequence(dim(Q)[1])
	sources.original <- sources.original[-to]
	sources.new <- sources.original
	sources.new[which(sources.new>to)] <- sources.new[which(sources.new>to)] + num.sources - 1
	for (i in sequence(num.sources)) {
		Q.expanded[sources.new[i], to + i - 1] <- Q[sources.original[i], to]
	}
	starting <- rep(0, dim(Q.expanded)[1])
	starting[from.new] <- 1
	diag(Q.expanded) <- 0
	diag(Q.expanded) <- -rowSums(Q.expanded)
	
	for (i in sequence(num.sources)) {
		colnames(Q.expanded)[to + i - 1] <- paste(colnames(Q)[to], "_from_", colnames(Q)[sources.original[i]], sep="")
	}

	
	P <- expm(t(Q.expanded) * total.time, method=c("Ward77")) %*% starting
	P <- P[to:(to+num.sources-1),]
	return(P)
}

PlotBubbleMatrix <- function(Q, main="", special=Inf, cex=1){
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

PlotTransitionNetwork <- function(Q, main="", layout.fn = layout.circle, ...) {
	require(igraph)
	diag(Q) <- 0
	Q<-Q/max(Q)
	g <- graph.adjacency(Q, weighted=TRUE, mode="directed")
	g.layout <- layout.fn(g)
	plot(g, layout=g.layout, edge.width=10*get.edge.attribute(g, "weight"), edge.curved=TRUE, ...)
}




phen <- c(0.0511,0.0531,0.0108,0.0084,0.0116,0.0006,0.0116) 
names(phen)<-c("EEtoEU","DUtoEU","DEtoEE","EUtoDU","EEtoDE","EUtoDE","EEtoDE")
mm<-matrix(,6,6)
mm[2,1]<-1
mm[3,1]<-2
mm[1,2]<-3
mm[1,3]<-4
mm[2,4]<-5
mm[1,5]<-6
mm[3,6]<-7
rate<-mm
#Simple bookkeeping to ensure things go in the right place:
rate[is.na(rate)]=8
Q<-matrix(0,6,6)
Q[] <- c(phen, 0)[rate]
diag(Q) <- -rowSums(Q)
#Fixes the starting state (i.e., unity=1):        
liks<-c(1,0,0,0,0,0)
pathways<-(expm(t(Q) * 10000000, method=c("Ward77")) %*% liks)[4:6]
names(pathways)<-c("Climate", "Simul", "Trait")
print(pathways)

Q.new <- Q
Q.new[1,4] <- Q.new[1,5]
Q.new[3,4] <- Q.new[3,6]
Q.new[-5,]->Q.new
Q.new[-5,]->Q.new
Q.new[,-5] ->Q.new
Q.new[,-5] ->Q.new
rownames(Q.new) <- c("A", "B", "C", "D")
colnames(Q.new) <- c("A", "B", "C", "D")
ComputePathwayProbability(Q.new, 1, 4)