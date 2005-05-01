\name{fwdmvPairsPlot}
\alias{fwdmvPairsPlot}

\title{A Pairs-like Plots of an fwdmv Object}

\description{
For initial \code{fwdmv} objects this function produces a plot almost identical to the function \code{pairs}.  If tentative groups have been assigned (either through the use of the \code{partition} function or directly in the call to \code{fwdmv}) then the units in each tentative group are plotted with different symbols.  Additionally, the indicies of the unassigned units are plotted.
}

\usage{
fwdmvPairsPlot(x)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
}

\value{
an empty list is returned invisibly.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)

g1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 53, 55, 56)

g2 <- c(57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)

fondi.fwdmv <- fwdmv(fondi.dat, groups = list(g1, g2))

fwdmvPairsPlot(fondi.fwdmv)
}

\keyword{hplot}
