\name{plot.fwdtr.test}
\alias{plot.fwdtr.test}

\title{ Plot method for fwdtr.test }

\description{
This function is the generic plot method for fwdtr.test objects.  It produces a series of fan plots for a set of variables.  It allows the user to select for which variables to view the fan plot.
}

\usage{
\method{plot}{fwdtr.test}(x, psfrag.labels = FALSE, ...)
}

\arguments{
  \item{x}{a \code{fwdtr.test} object.}
  
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}


  \item{\dots}{these arguments will be passed on to the selected plot.}
}

\details{
 It is also possible to have the confirmatory fan plots for each variable as panels in a single page
}

\value{
 the value returned by the selected plot.
}

\references{ Atkinson, A. C., Riani, M. and Cerioli, A. (2004)  Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York.}

\author{ Fabrizio Laurini }

\seealso{\code{\link{fwdtr.object}}, \code{\link{profile.fwdtr.object}}}

\keyword{hplot}
