\name{fwdtrProfilePlot}
\alias{fwdtrProfilePlot}

\title{Plot method for profile.fwdtr objects}
\description{
Produces a plot of profile loglikelihoods of transformation parameters in a particular step of the forward search
}
\usage{
fwdtrProfilePlot(x, psfrag.labels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{profile.fwdtr} object. }
  \item{psfrag.labels}{a logical value.  If \code{TRUE} then the x, y, and main labels are set to "xlab", "ylab", and "main" for replacement via the psfrag utility.}
}
\details{
This plot shows which variables show a sharp definition of the estimates of the transformation parameters and which are the variables whose value of lambda is not very well determined.}
\value{
 an empty list is invisibly returned.
}
\references{ Atkinson, A. C., Riani, M. and Cerioli, A. (2004)  Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York.}
\author{Fabrizio Laurini}


\seealso{\code{\link{fwdtr.object}}}
\examples{data(mussels.dat)

## Forward search on untransformed data 
## Compute max lik. estimates of tramsformation parameters

l.mle <- fwdtr(mussels.dat)

## Compute profile loglikelihoods for transformation parameters in the
## last step of the search and create an object of class profile.fwdtr

l.profile.mle <- profile.fwdtr(l.mle)

## plot the profile loglikelihoods of transformation parameters for
## each variable 

fwdtrProfilePlot(l.profile.mle)

## Test variables 2 and 5
## The forward is based on untransformed data for variables 1, 3 and 4
## sqrt for variable 2 and third root for variable 5

l.mle <- fwdtr(mussels.dat, lambda = c(0.5, 1/3), col.to.transform = c(2, 5))

## build profile likelihood for transformed variables.
## Profile function takes into account that we have also untransformed variables.

l.profile.mle <- profile.fwdtr(l.mle)
fwdtrProfilePlot(l.profile.mle)
}

\keyword{hplot}
