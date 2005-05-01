\name{fwdmvPrincompPlot}
\alias{fwdmvPrincompPlot}

\title{Plot the Principal Components from an fwdmv Object}

\description{
Plots (for each group in the fwdmv object) the principal components against the subset size.
}

\usage{
fwdmvPrincompPlot(x, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
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
fwdmvPrincompPlot(fondi.1)
}

\keyword{hplot}
