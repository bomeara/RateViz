\name{PlotBubbleMatrix}
\alias{PlotBubbleMatrix}
\title{
PlotBubbleMatrix
}
\description{
Plots a rate matrix such that circle area correlates with the rate
}
\usage{
PlotBubbleMatrix(Q, main = "", special = Inf, cex = 1)
}
\arguments{
  \item{Q}{
Instantaneous rate matrix
}
  \item{main}{
Optional plot title
}
  \item{special}{
If you want to highlight a particular element, input its number here
}
  \item{cex}{
Adjustment for text size
}
}
\value{
None
}

\author{
Brian C. O'Meara
}

\examples{
#A random sample transition matrix
Q <- matrix(runif(9), 3, 3)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
colnames(Q) <- c("A", "B", "C")
rownames(Q) <- c("A", "B", "C")
print(Q)
PlotBubbleMatrix(Q)
}