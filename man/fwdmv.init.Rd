\name{fwdmv.init}
\alias{fwdmv.init}
\alias{fwdmv.init.pureR}

\title{Multivariate Forward Search for Ungrouped Data}

\description{
This function computes a multivariate forward search for ungrouped data.  Several diagnostic statistics are monitored during the search: see \code{fwdmv.object}.  Note that this function is called by \code{\link{fwdmv}} when no tentative groups are specified.  It is recommended that \code{\link{fwdmv}} be used for all multivariate forward searches.
}

\usage{
fwdmv.init(X, bsb = ellipse.subset, n.bsb, scaled = TRUE, monitor = "all")
}

\arguments{
  \item{X}{a matrix or data frame containing the multivariate data set.}
  \item{bsb}{a function of two variables: the multivariate data in matrix form \code{X} and the number of units in the initial subset \code{size}.  The default \code{bsb = ellipse.subset} computes the initial subset using robustly centered ellipses.  Other choices include \code{bsb = bb.subset} to compute the initial subset using bivariate boxplots, \code{bsb = mcd.subset} to compute the initial subset using mcd distances, and \code{bsb = random.subset} for a randomly determined initial subset.  Alternatively, the initial subset my be specified directly by providing an integer vector containing the indices of the units to be in the initial subset.}
  \item{n.bsb}{a positive integer specifying the size of the initial subset; ignored if \code{bsb} is not a function.}
  \item{scaled}{a logical value.  If \code{TRUE} then scaled Mahalanobis distances are used during the forward search.}
  \item{monitor}{a character vector specifying which statistics are to be monitored during the forward search.  The default value "all" monitors all statistics.  Otherwise choose from "distance", "center", "cov", "determinant", "unit", "max", "mth", "min", and "mpo".}
}

\details{
This function computes the Forward Search as described in chapter 3 of ARC.  The initial subset can be specified directly in the argument \code{bsb} or computed from the data.  By default \code{bsb} is a function for computing the initial subset using robustly centered ellipses.  Given a subset of \code{m} units the next subset is the \code{m+1} units with smallest Mahalanobis distances calculated using the center and covariance matrix estimates computed from the units currently in the subset.  This process is repeated until the subset contains all of the units and diagnostic statistics are monitored for each subset.
}

\value{
  an \code{fwdmv} object.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv}}, \code{\link{fwdmv.object}}}

\examples{
data(fondi.dat)
fondi.init <- fwdmv.init(fondi.dat)

## find the intial subset using robust bivariate ellipses ##
## start with an initial subset size of 17 units          ##
fondi.init <- fwdmv.init(fondi.dat,bsb=ellipse.subset(fondi.dat,17))

## find the intial subset using bivariate boxplots and ##
## start with an initial subset size of 17 units       ##
fondi.init <- fwdmv.init(fondi.dat,bsb=bb.subset(fondi.dat,17))

}

\keyword{multivariate}
