\name{mcd.subset}
\alias{mcd.subset}

\title{Initial Subset by MCD Distances}

\description{
This function uses the robust distances obtained from \code{cov.mcd} to compute the initial subset for the multivariate forward search.  It is not intended that users should call this function directly.}

\usage{
mcd.subset(X, size)
}

\arguments{
  \item{X}{a numeric matrix containing the multivariate data set.}
  \item{size}{an integer specifying the size of the initial subset.}
}

\details{
The function \code{cov.mcd} in package \code{lqs} is used to robustly estimate the center and covariance matrix of \code{X}.  Robust Mahalanobis distances are computed and the intial subset is taken to be the \code{size} units with the smallest robust distances.
}

\value{
an integer vector with length \code{size} containing the initial subset.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\examples{
data(fondi.dat)

fondi.fwdmv <- fwdmv(fondi.dat, bsb = mcd.subset)
}

\keyword{utilities}
