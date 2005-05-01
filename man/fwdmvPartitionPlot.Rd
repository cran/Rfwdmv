\name{fwdmvPartitionPlot}
\alias{fwdmvPartitionPlot}

\title{Plot a Partitioned fwdmv Object}

\description{
Produces a plot of the Mahalanobis distances similar to that produced by \code{fwdmvDistancePlot}.  The trajectories for units assigned to tentative groups are not drawn.  Instead, the median distance is drawn for each group in the \code{fwdmv} object.}

\usage{
fwdmvPartitionPlot(x, pts = NULL, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{pts}{optionally include \code{x$pts} to draw the segment selected in \code{partition}.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\value{
an empty list is invisibly returned.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{partition}}, \code{\link{fwdmv}}, \code{\link{fwdmv.object}}}

\examples{
# data(fondi.dat)
# fondi.init <- fwdmv.init(fondi.dat)
# p1 <- partition(fondi.init)

# draw a line segment intersecting several trajectories #
# view the allocation #

# fwdmvPartitionPlot(p1)
}

\keyword{hplot}
