\name{fwdmvEccentricityPlot}
\alias{fwdmvEccentricityPlot}

\title{Eccentricity Plot}

\description{
Plots (for each group in the \code{fwdmv} object) the fraction of variance explained by the eigenvalues in \code{which} and the eccentricity of the ellipse from the same two eigenvalues.
}

\usage{
fwdmvEccentricityPlot(x, which = c(1, 2), psfrag.labels = FALSE)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{which}{an integer vector of length 2.  Specifies the eignevalues to use in the plot.}
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
fwdmvEccentricityPlot(fondi.1)
}

\keyword{hplot}
