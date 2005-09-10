\name{identify.fwdmvRandomStart}
\alias{identify.fwdmvRandomStart}

\title{Identify Starting Subsets From a Random Starts Forward Plot}

\description{
Allows the user to select a trajectory from the Random Starts forward plot, highlights it, and reports the starting subset(s) that lead to the trajectory.
}

\usage{
\method{identify}{fwdmvRandomStart}(x, ...)
}

\arguments{
  \item{x}{an \code{fwdmvRandomStart} object created by \code{\link{fwdmvRandomStart}}.}
  \item{\dots}{additional argument required by the geneic \code{\link{identify}} function.}
}

\value{
  a numeric matrix whose rows contain the starting subsets that lead to the user selected trajectory.
}

\seealso{\code{\link{fwdmvRandomStart}}}

\examples{
data(bank.dat)
bank.rfs <- fwdmvRandomStart(bank.dat)
# use identify and plot on bank.rfs #
}

\keyword{iplot}

