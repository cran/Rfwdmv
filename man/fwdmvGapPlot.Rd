\name{fwdmvGapPlot}
\alias{fwdmvGapPlot}

\title{Plot the Gap in an fwdmv Object}

\description{
Plots the minimum Mahalanobis distance among points not in the subset minus the maximum Mahalanobis distance among points in the subset the (m+1)th ordered Mahalanobis distance minus the mth ordered distance.
}

\usage{
fwdmvGapPlot(x, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\value{
  an empty list is returned invisibly.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fwdmvGapPlot(fondi.1)
}

\keyword{hplot}
