\name{partition}
\alias{partition}

\title{Interactive Group Assignment}

\description{
This function produces a plot of the Mahalanobis distances versus subset size similar to that produced by \code{fwdmvDistancePlot}.  The function then waits for the user to draw a line segment on the plot (which is done with two consecutive mouse clicks).  The units whose trajectories cross this line segment are assigned to the group specified by the \code{group} argument.
}

\usage{
partition(x, group = "next") 
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{group}{an integer value specifying which group the units corresponding to the selected trajectories should be assigned.  Possible values are \code{1} to \code{n.groups + 1} where \code{n.groups} is the number of groups in the \code{fwdmv} object.  The default \code{"next"} is to assign units to a new group.}
}

\value{
  an \code{fwdmv} object with the \code{groups} element and the \code{unassigned} element updated according to the description.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}, \code{\link{fwdmvPartitionPlot}}}

\examples{
# data(fondi.dat)
# fondi.init <- fwdmv.init(fondi.dat)
# p1 <- partition(fondi.init)

# draw a line segment intersecting several trajectories #
# and view the allocation #
}

\keyword{hplot}
\keyword{iplot}
