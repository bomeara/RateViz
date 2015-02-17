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
%%     ~~Describe \code{Q} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{correct.diag}{
%%     ~~Describe \code{correct.diag} here~~
}
  \item{total.time}{
%%     ~~Describe \code{total.time} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
ComputePathwayProbability(Q)
