\name{fwdmvEntryPlot}
\alias{fwdmvEntryPlot}

\title{Plot the Entry Order in an fwdmv Object}

\description{
Produces an entry order plot of an fwdmv object.
}

\usage{
fwdmvEntryPlot(x, entry.order = "first", subset.size = -1, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{entry.order}{a character vector of length 1 specifying how the rows should be ordered.  The possibilities are "first", "final", "natural" and "integer".  If \code{entry.order == "integer"} then \code{entry.order.n} must be supplied as well.}
  \item{subset.size}{an integer value giving the subset size to be used when \code{entry.order == "integer"}.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}

\details{
An entry oder plot is a matrix where the (i,j) cell is black if the unit represented by row i is in the subset during step j and white otherwise.  The rows can be ordered in four ways: (1) First entry order: the rows are ordered (from bottom to top) by the first time the unit enters the subset; (2) Final entry order: the row are ordered (from bottom to top) by the last time the unit enters the subset; (3) Natural entry order: the rows appear (from bottom to top) in the same order as in the data; and (4) Integer entry order: the rows are ordered (from bottom to top) in the same order as the distances in the subset of \code{subset.size} units.
}

\value{
  an empty list is returned invisibly.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fwdmvEntryPlot(fondi.1)
}

\keyword{hplot}
