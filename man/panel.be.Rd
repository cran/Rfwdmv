\name{panel.be}
\alias{panel.be}

\title{Bivariate Ellipse Panel Function}

\description{
This function is intended to be used in the panel argument of \code{fwdmvPrePlot}.  It computes a bivariate ellipse for the given panel.
}

\usage{
panel.be(x, y, scale = 1)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric vector with the same length as \code{x}.}
  \item{scale}{a positive numeric value to scale the ellipse.}
}

\value{
a list with elements x and y giving the points of the ellipse.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmvPrePlot}}}

\examples{
data(heads.dat)

fwdmvPrePlot(heads.dat, panel = panel.be)
}

\keyword{hplot}

