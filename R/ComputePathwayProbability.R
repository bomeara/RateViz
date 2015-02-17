ComputePathwayProbability <-
function(Q, from, to, correct.diag=TRUE, total.time=1000000) {
	
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
