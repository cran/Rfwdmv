\name{fwdmv}
\alias{fwdmv}

\title{Multivariate Forward Search}

\description{
This function computes a multivariate forward search.  Several diagnostic statistics are monitored during the search: see \code{\link{fwdmv.object}}.
}

\usage{
fwdmv(X, groups = NULL, alpha = 0.6, beta = 0.75, bsb = ellipse.subset, balanced = TRUE, scaled = TRUE, constrained = TRUE, monitor = "all")
}

\arguments{
  \item{X}{a matrix or data frame containing a multivariate data set.}
  \item{groups}{a list of one or more integer vectors specifying the tentative groups.  All elements must be unique.  Units not belonging to any group are classified as unassigned.  If omitted, all of the data are assumed to come from a single multivariate normal population.}
  \item{alpha}{a numeric value between 0 and 1 specifying the fraction of the units in each group that will be included in the initial subset.}
  \item{beta}{a numeric value between \code{alpha} and 1 specifying the fraction of the units in each tentative group that must be included in the subset before the unassigned units are allowed to be included.  A large value of \code{beta} insures that the centroid and variance-covariance matrix estimates stabilize before the unassigned units enter the subsets.}
  \item{bsb}{a function of two variables: the multivariate data in matrix form \code{X} and the number of units in the initial subset \code{size}.  The default \code{bsb = ellipse.subset} computes the initial subset using robustly centered ellipses.  Other choices include \code{bsb = bb.subset} to compute the initial subset using bivariate boxplots, \code{bsb = mcd.subset} to compute the initial subset using mcd distances, and \code{bsb = random.subset} for a randomly determined initial subset.  Alternatively, the initial subset my be specified directly by providing an integer vector containing the indices of the units to be in the initial subset.}
  \item{balanced}{a logical value.  If \code{TRUE} then units are added to the subset so that the group ratios in the subset stay as close as possible to the group ratio in the data.}
  \item{scaled}{a logical value.  If \code{TRUE} then the Mahalanobis distances are scaled using the 2p root of the determinant.  This is intended to compensate for clusters with significantly different dispersions.}
  \item{constrained}{a logical value.  If TRUE then the forward search chooses units from the tentative groups until each group is \code{beta} full; then unassigned units are allowed into the subset.  If FALSE then unassigned units may enter the subset at anytime during the forward search.  Note that when \code{constrained == F}, the argument \code{balanced} is ignored.}
  \item{monitor}{a character vector specifying which statistics are to be monitored during the forward search.  The default value "all" monitors all statistics.  Otherwise choose from "distance", "center", "cov", "determinant", "unit", "max", "mth", "min", "mpo", "nearest", and "misclassified".}
}

\details{
Initial group subsets of size \code{alpha * nbsb[i]} (where \code{nbsb[i]} is the number of units assigned to tentative group i) are obtained by running the initialization function on each group.  Estimates of the center and covariance matrix are computed for each group using the units currently in the group subset.  The Mahalanobis distance for each unit in a tentative group is computed using the center and covariance matrix estimates for that group.  The Mahalanobis distance for each unassigned unit is computed by calculating the distance to each group and taking the minimum.  If the search is balanced then one unit is added to the subset that is currently the furthest below the population ratio.  If the search is not balanced then the unit (not in any subset) with the smallest distance is allocated to the nearest group.  If the search is constrained then the unassigned units are not allowed into the group subsets until each group subset contains a fraction \code{beta} of the units in the tentative groups.  If the search is not constrained then the unassigned units may enter the subset at any time during the search.
}

\value{
  a list with class \code{fwdmv}.
}

\references{
Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}}

\examples{
data(fondi.dat)

g1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 53, 55, 56)

g2 <- c(57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)

fondi.fwdmv <- fwdmv(fondi.dat, groups = list(g1, g2))
}

\keyword{multivariate}
