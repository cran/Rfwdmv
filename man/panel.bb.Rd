\name{panel.bb}
\alias{panel.bb}

\title{Bivariate Box Plots Panel Function}

\description{
This function is intended to be used in the panel argument of \code{fwdmvPrePlot}.  It computes a bivariate boxplot for the given panel.
}

\usage{
panel.bb(x, y, scale = 1)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric vector with the same length as \code{x}.}
  \item{scale}{a positive numeric value to scale the boxplot contour.}
}

\value{
a list with elements x and y giving the points of the computed contour.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmvPrePlot}}}

\examples{
data(heads.dat)

fwdmvPrePlot(heads.dat, panel = panel.bb)
}

\keyword{hplot}

