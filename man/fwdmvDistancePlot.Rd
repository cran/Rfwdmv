\name{fwdmvDistancePlot}
\alias{fwdmvDistancePlot}

\title{Plot the Mahalanobis Distances in an fwdmv Object}

\description{
Plots the Mahalanobis distance against the subset size for each unit in the \code{fwdmv} object.
}

\usage{
fwdmvDistancePlot(x, group = NULL, id = FALSE, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
	\item{group}{show only the trajectories for the units in this group.}
  \item{id}{identify trajectories in the plot.  If \code{id = TRUE} then trajectories in the plot can be identified interactively with the mouse.  If \code{id} is an integer vector then the trajectories for those units are identified in the plot.  If \code{id = "all"} then all of the trajectories are labelled.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\value{
  a list is invisibly returned.  If \code{id = TRUE} then the list has an element named \code{selected} containing the indices of the identified trajectories.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)

fwdmvDistancePlot(fondi.1)

# Use 'id = TRUE' for interactive trajectory identification.

fwdmvDistancePlot(fondi.1, id = c(39, 52, 96))

g1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 53, 55, 56)

g2 <- c(57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)

fondi.2 <- fwdmv(fondi.dat, groups = list(g1, g2))

fwdmvDistancePlot(fondi.2, group = 2)
}

\keyword{hplot}
\keyword{iplot}
