\name{fwdtrMlePlot}
\alias{fwdtrMlePlot}

\title{ Plot method for fwdtr objects }
\description{
Produces a plot of maximum likelihood estimates of transformation parameters during all steps of the forward search.
}
\usage{
fwdtrMlePlot(x, psfrag.labels = FALSE)
}

\arguments{
  \item{x}{ a \code{fwdtr} object.}
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}
\details{
Estimates that vary wildly are associated with variables that do no have to be transformed.
}
\value{
an empty list is invisibly returned
}
\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}
\author{Fabrizio Laurini}

\seealso{\code{\link{fwdtr.object}}}
\examples{
## Forward search on untrasformed data 

 data(mussels.dat)
 l.mle<-fwdtr(mussels.dat)

## Plot maximum likelihood estimates of the trasnformation parameters
## in each step of the forward search 

fwdtrMlePlot(l.mle)

 
## Test variables 2 and 5
## The forward is based on untransformed data for variables 1, 3 and 4
## sqrt for variable 2 and third root for variable 5

 l.mle <- fwdtr(mussels.dat, lambda = c(0.5, 1/3), col.to.transform = c(2,5))

## plot trajectories of maximum likelihood estimates of transformation parameters
## for variables 2 and 5

fwdtrMlePlot(l.mle)
}

\keyword{hplot}
