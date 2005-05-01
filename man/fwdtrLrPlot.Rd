\name{fwdtrLrPlot}
\alias{fwdtrLrPlot}

\title{ Plot method for fwdtr objects }
\description{
   Produces a plot of the likelihood ratio of transformation parameters during all steps of the forward search.
}
\usage{
fwdtrLrPlot(x, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{a \code{fwdtr} object. }
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}
\details{
The horizontal lines drawn in the plot refer to the 95\% and 99\% quantiles of the associated chisquare distribution.
}
\value{
an empty list is invisibly returned.
}
\references{
Atkinson, A. C., Riani, M. and Cerioli, A. (2004)  Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York.
}
\author{Fabrizio Laurini}


\seealso{\code{\link{fwdtr.object}}, \code{\link{profile.fwdtr.object}}}
\examples{
## Forward search on untrasformed data 

 data(mussels.dat)
 l.mle <- fwdtr(mussels.dat)
 
## plot the likelihood ratio test 

 fwdtrLrPlot(l.mle)
}

\keyword{hplot}
