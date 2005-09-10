\name{fwdmvRandomStart}
\alias{fwdmvRandomStart}

\title{Random Starting Points for the Multivariate Forward Search}

\description{
Produces a forward plot of the trajectories of the minimum distance outside the subset for several multivariate forward searches.  Peaked trajectories are often indicative of clusters in the data.  The identify function can be used to highlight these trajectories in the plot and to obtain one or more starting subsets that are quite likely to reveal group structure in the data.
}

\usage{
fwdmvRandomStart(X, n.starts = 50, scaled = TRUE, plot.it = TRUE)
}

\arguments{
  \item{X}{a matrix or data frame containing the data.}
  \item{n.starts}{a positive integer specifying the number of random starting subsets.}
  \item{scaled}{a logical flag.  If \code{TRUE} scaled mahalanobis distances are used.}
  \item{plot.it}{a logical flag.  If \code{TRUE} the minimum distance outside the subset trajectories are plotted.}
}

\details{
It is necessary to store the object returned by this function in order to use the identify function.  Additionally, the plot method can be used to redraw the plot.
}

\value{
  A list with class \code{fwdmvRandomStart} containing the following components:

  \item{starts}{a numeric matrix whose rows contain the randomly selected starting subsets.}
  \item{mpomat}{a numeric matrix where each column contains the minimum distance outside the subset trajectory for a particular multivariate forward search.}
  \item{plot.domain}{a numeric vector containing the subset sizes used during the forward search.}
}

\seealso{
  \code{\link{fwdmv}},
  \code{\link{identify.fwdmvRandomStart}},
  \code{\link{plot.fwdmvRandomStart}}
}

\examples{
data(bank.dat)
bank.rfs <- fwdmvRandomStart(bank.dat)
# use identify and plot on bank.rfs #
}

\keyword{multivariate}

