\name{PlotTransitionNetwork}
\alias{PlotTransitionNetwork}
\title{
PlotTransitionNetwork
}
\description{
This plots a network showing transition rates between given states. It uses igraph for plotting.
}
\usage{
PlotTransitionNetwork(Q, main = "", layout.fn = layout.circle, ...)
}
\\arguments{
  \item{Q}{
Instantaneous rate matrix
}
  \item{main}{
Title for plot (optional)
}
  \item{layout.fn}{
Layout function from igraph. By default this uses a circle layout, but you may use other igraph functions. layout.fruchterman.reingold may generate a good looking graph.}
  \item{\dots}{
Other arguments to pass to plot.igraph()}
}

\value{
Nothing is returned. Except beauty.
}
\author{
Brian C. O'Meara
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#A random sample transition matrix
Q <- matrix(runif(9), 3, 3)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
colnames(Q) <- c("A", "B", "C")
rownames(Q) <- c("A", "B", "C")
print(Q)
PlotTransitionNetwork(Q)
