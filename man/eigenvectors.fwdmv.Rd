\name{eigenvectors.fwdmv}
\alias{eigenvectors.fwdmv}

\title{Retrieve an Eigenvector of the Covariance Matrix Estimate from an fwdmv Object}

\description{
An accessor method to retrieve an eigenvector of the covariance matrix estimate from an \code{fwdmv} object.
}

\usage{
eigenvectors.fwdmv(x, which.vector = 1)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{which.vector}{an integer specifying which eigenvector should be returned.}
}

\value{
a list with one element for each group in the \code{fwdmv} object \code{x}.  Each element is a matrix where row i contains the components of the \code{which.vector} eigenvector of the covariance matrix estimate computed during step i of the forward search.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fondi.evec1 <- eigenvectors.fwdmv(fondi.1, which.vector = 1)
fondi.evec3 <- eigenvectors.fwdmv(fondi.1, which.vector = 3)
}

\keyword{methods}
