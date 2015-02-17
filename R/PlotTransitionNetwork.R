PlotTransitionNetwork <-
function(Q, main="", layout.fn = layout.circle, ...) {
	diag(Q) <- 0
	Q<-Q/max(Q)
	g <- graph.adjacency(Q, weighted=TRUE, mode="directed")
	g.layout <- layout.fn(g)
	plot(g, layout=g.layout, edge.width=10*get.edge.attribute(g, "weight"), edge.curved=TRUE, ...)
}
