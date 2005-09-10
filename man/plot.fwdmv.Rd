\name{plot.fwdmv}

\alias{plot.fwdmv}

\title{Plot method for fwdmv objects.}

\description{
This function is the generic plot method for \code{fwdmv} objects.  It allows the user to select from a subset of the plotting functions in the Rfwdmv package.
}

\usage{
\method{plot}{fwdmv}(x, ...)
}

\arguments{
  \item{x}{an fwdmv object.}
  \item{\dots}{these arguments will be passed on to the selected plot function.}
}

\details{
Plot functions that require additional arguments must be called directly and are thus not available from this function: see \code{\link{fwdmv.object}} for a list of all available plot methods.
}

\value{
  the value returned by the selected plot.
}

\author{Kjell Konis}

\seealso{ \code{\link{fwdmv}},
          \code{\link{fwdmvChangePlot}},
          \code{\link{fwdmvCovariancePlot}},
          \code{\link{fwdmvDeterminantPlot}},
          \code{\link{fwdmvDistancePlot}},
          \code{\link{fwdmvPrincompPlot}},
          \code{\link{fwdmvEntryPlot}},
          \code{\link{fwdmvGapPlot}},
          \code{\link{fwdmvMinmaxPlot}},
          \code{\link{fwdmvPairsPlot}} }

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
plot(fondi.1)
}

\keyword{hplot}
