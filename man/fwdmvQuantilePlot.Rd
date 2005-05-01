\name{fwdmvQuantilePlot}
\alias{fwdmvQuantilePlot}

\title{Plot Trajectories over Quantiles of the Distances}

\description{
Produces a 6 panel plot.  The first panel contains the Mahalanobis distances for the units in the subset of size \code{subset.size}.  The backgrounds of the 5 remaining plots are quantiles of these distances.  The dark lines in the 5 remaining plots are the trajectories of the units entering the subset during each of the next 5 steps of the forward search.
}

\usage{
fwdmvQuantilePlot(x, subset.size, probs = "default", page = 1)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{subset.size}{an integer specifying (by its size) the subset to be used in the first panel.}
  \item{probs}{an ordered numeric vector of probabilities used by the quantile function.}
  \item{page}{an integer specifying the page.  When \code{page == 1} the trajectories for the next 5 steps are drawn, when \code{page == 2} the next 5 steps are skipped and trajectories for the 5 following step are drawn, etc.}
}

\details{
when \code{probs == "default"} the following quantiles are computed: 0.025, 0.05, 0.125, 0.25, 0.50, 0.75, 0.875, 0.95, 0.975.
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
fwdmvQuantilePlot(fondi.1, subset.size = 45)
}

\keyword{hplot}
