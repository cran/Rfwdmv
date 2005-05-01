\name{fwdmvEigenvectorPlot}
\alias{fwdmvEigenvectorPlot}

\title{Plot the Components of an Eigenvector in an fwdmv Object}

\description{
Plots (for each group in the \code{fwdmv} object) the components of the specified eigenvector of the covariance matrix against the subset size.
}

\usage{
fwdmvEigenvectorPlot(x, which.vector = 1, correlation = FALSE, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{which.vector}{an integer value used to select the eigenvector whose components are to be plotted.}
	\item{correlation}{a logical value.  If \code{TRUE} then the eigenvectors of the correlation matrix are plotted.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
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
fwdmvEigenvectorPlot(fondi.1)
}

\keyword{hplot}
