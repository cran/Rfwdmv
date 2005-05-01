\name{fwdmvPrePlot}
\alias{fwdmvPrePlot}

\title{Pairs Plot in Rfwdmv}

\description{
This function produces a pairs plot with sumperimposed contours.
}

\usage{
fwdmvPrePlot(X, panel = panel.be, plot.diagonal = TRUE)
}

\arguments{
  \item{X}{a matrix or data frame.}
  \item{panel}{a function for computing the contours.  The Rfwdmv package includes \code{panel.be} for bivariate ellipses, \code{panel.me} for median centered bivariate ellipses, and \code{panel.bb} for bivariate box plots.}
  \item{plot.diagonal}{a logical value.  If \code{TRUE} then univariate boxplots are drawn along the diagonal.}
}

\value{
an empty list is invisibly returned.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\note{ The bivariate boxplots calculated from B-splines provide a useful tool for a preliminary examination of the data. 
The non elliptical shape of the countours is an indication or non normality.}

\seealso{\code{\link{panel.be}}, \code{\link{panel.me}}, \code{\link{panel.bb}}}

\examples{
data(heads.dat)

fwdmvPrePlot(heads.dat, panel = panel.be)
fwdmvPrePlot(heads.dat, panel = panel.me)
fwdmvPrePlot(heads.dat, panel = panel.bb)
}

\keyword{hplot}

