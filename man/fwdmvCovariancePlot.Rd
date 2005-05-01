\name{fwdmvCovariancePlot}
\alias{fwdmvCovariancePlot}

\title{Plot the Covariance Matrix of an fwdmv Object}

\description{
Plots the elements of the covariance matrix of each group against the subset size.
}

\usage{
fwdmvCovariancePlot(x, id = FALSE, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{id}{a logical value.  If \code{TRUE} then the curves in the plot can be identified interactively with the mouse.  Alternatively, \code{id = "all"} will identify all of the curves in the plot.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\value{
  a list is invisibly returned.  If \code{id = TRUE} then the list has an element named \code{selected} containing the identified curves.  The curves are identified by their column index in the fwdmv object.  If there is more than one group then the numbers 1 through p are for the first group, p+1 through 2p for the second group and so on.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)

g1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 53, 55, 56)

g2 <- c(57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)

fondi.fwdmv <- fwdmv(fondi.dat, groups = list(g1, g2))

fwdmvCovariancePlot(fondi.fwdmv)

fwdmvCovariancePlot(fondi.fwdmv, id = "all")

# Use 'id = TRUE' for interactive curve labels.
}

\keyword{hplot}
\keyword{iplot}
