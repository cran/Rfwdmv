\name{fwdmvMinmaxPlot}
\alias{fwdmvMinmaxPlot}

\title{Plot Minimum and Maximum Distances in an fwdmv Object}

\description{
A two panel plot.  The first panel shows the maximum Mahalanobis distance among units in the subset and mth ordered Mahalanobis distance.  The second panel shows the (m+1)th ordered Mahalanobis distance and the minimum Mahalanobis among units in the complement of the subset.
}

\usage{
fwdmvMinmaxPlot(x, psfrag.labels = FALSE)
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
fwdmvMinmaxPlot(fondi.1)
}

\keyword{hplot}
