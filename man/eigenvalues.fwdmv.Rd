\name{eigenvalues.fwdmv}
\alias{eigenvalues.fwdmv}

\title{Compute the Eigenvalues from an fwdmv Object}

\description{
An accessor method to retrieve the eigenvalues of the covariance matrix estimates from an \code{fwdmv} object.
}

\usage{
eigenvalues.fwdmv(x)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
}

\value{
a list with one element for each group in the \code{fwdmv} object \code{x}.  Each element is a matrix where row i contains the eigenvalues of the covariance matrix estimate computed during step i of the forward search.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fondi.evals <- eigenvalues.fwdmv(fondi.1)
}

\keyword{utilities}
