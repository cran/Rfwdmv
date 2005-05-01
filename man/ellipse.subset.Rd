\name{ellipse.subset}
\alias{ellipse.subset}

\title{Initial Subset by Robustly Centered Ellipses}

\description{
Computes the initial subset using robustly centered ellipses as described in ARC.
}

\usage{
ellipse.subset(X, size)
}

\arguments{
  \item{X}{a numeric matrix containing the multivariate data set.}
  \item{size}{an integer specifying the size of the initial subset.}
}

\references{
Riani M. and Zani S. (1997), An Iterative Method for the Detection of 
Multivariate Outliers", Metron, pp. 101-117.

Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)

fondi.fwdmv <- fwdmv(fondi.dat, bsb = ellipse.subset)
}

\keyword{utilities}
