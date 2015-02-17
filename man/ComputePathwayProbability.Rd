\name{ComputePathwayProbability}
\alias{ComputePathwayProbability}
\title{
ComputePathwayProbability}
\description{
This implements the rate summary used in Zanne et al. 2014. Using the rate matrix, one can calculate }
\usage{
ComputePathwayProbability(Q, from, to, correct.diag = TRUE, total.time = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
Instantaneous rate matrix
}
  \item{from}{
Index (position) of start state in the matrix (i.e., the second character)
}
  \item{to}{
Index (position) of end state in the matrix
}
  \item{correct.diag}{
Normally the diagonal is the negative sum of the non-diagonal elements in a row. If TRUE, this verifies that this is so
}
  \item{total.time}{
A large enough time such that there is essentially no chance of remaining unchanged}
}
\value{
A vector with the probability of each path into the final state.
}
\references{
Zanne et al. 2014
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
ComputePathwayProbability(Q, from=1, to=3)
}