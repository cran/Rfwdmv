\name{fwdmvChangePlot}
\alias{fwdmvChangePlot}

\title{Change Plot of an fwdmv Object}

\description{
A matrix with one row for each unit in the data and one column for each step of the forward search.  If the distance of unit i increases during step j of the forward search then cell (i,j) is filled - otherwise it is empty.
}

\usage{
fwdmvChangePlot(x, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\value{
  an empty list is invisibly returned.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fwdmvChangePlot(fondi.1)
}

\keyword{hplot}

