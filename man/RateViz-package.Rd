\name{RateViz-package}
\alias{RateViz-package}
\alias{RateViz}
\docType{package}
\title{
Visualizations for instantaneous rate matrices
}
\description{
This package uses an instantaneous rate matrix (such as is used for models of DNA evolution, morphological character change over evolutionary history, and other areas) and summarizes it in various ways. These could include bubble plots, a network showing transitions, or calculations to figure out likely paths between two states.
}
\details{
\tabular{ll}{
Package: \tab RateViz\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-17\cr
License: \tab GPL (>=2)\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Brian C. O'Meara and Jeremy M. Beaulieu

Maintainer: Brian C. O'Meara <omeara.brian@gmail.com>
}


\examples{
#A random sample transition matrix
Q <- matrix(runif(9), 3, 3)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
colnames(Q) <- c("A", "B", "C")
rownames(Q) <- c("A", "B", "C")
print(Q)
ComputePathwayProbability(Q)
PlotTransitionNetwork(Q)

}
