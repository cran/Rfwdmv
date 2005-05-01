\name{profile.fwdtr}
\alias{profile.fwdtr}

\title{Profile log-likelihood estimates of transformation parameters}
\description{
 This function computes the profile loglikelihoods of transformation parameters.
}
\usage{
profile.fwdtr(fitted, step.fwd = NULL, conf = 0.95, bounds = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{  an object of class \code{fwdtr}.  This object can be created by function \code{fwdtr}.}
  \item{step.fwd}{an integer value.  The step of the forward search for which profile loglikelihoods must be created.  If \code{step.fwd = NULL} the profile loglikelihood is computed for the last step of the forward search.}
  \item{conf}{ scalar between 0 and 1 which defines the marginal confidence interval of lambda for each variable.  The default is 95\% confidence interval.}
  \item{bounds}{  a \code{2 x 1} numeric vector containing the lower and upper limit of the x axis for each profile loglikelihood.  The default is the value of \code{boundaries} in \code{fwdtr.object}.}
  \item{\dots}{Further parameters for profile method}
}
\details{
In order to compute the profile log-likelihoods the values of the parameters which are not being varied are kept at their maximum likelihood estimates. The loglikelihoods are roughly parabolic close to zero although not necessary log concave further away from the maximum. The confidence interval for each value of lambda is based on the asymptotic chi-square distribution of twice the loglikelihood ratio.
}
\value{
an object of class \code{profile.fwdtr}.
}
\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004)  Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York. }
\author{Fabrizio Laurini}


\seealso{\code{\link{profile.fwdtr.object}}, \code{\link{fwdtrProfilePlot}}  }
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

## Forward search on transformed data as specified in vector lamdba 

data(mussels.dat)
l.mle <- fwdtr(mussels.dat, lambda = c(1, 0.5, 1, 0, 1/3))

## Compute profile loglikelihoods for transformation parameters in the
## last step of the search and create an object of class profile.fwdtr

l.profile.mle <- profile.fwdtr(l.mle) 

## plot the profile loglikelihoods of transformation parameters for
## each variable

fwdtrProfilePlot(l.profile.mle)

## estimate a common value of lambda for all the variables and use 1/3
## to order Mahalanobis distances in each step of the search

data(mussels.dat)
l.mle <- fwdtr(mussels.dat, lambda = 1/3, one.lambda = TRUE)
l.profile.mle <- profile.fwdtr(l.mle)
fwdtrProfilePlot(l.profile.mle)

## Test variables 2 and 5
## The forward is based on untransformed data for variables 1, 3 and 4
## sqrt for variable 2 and third root for variable 5

l.mle <- fwdtr(mussels.dat, lambda = c(0.5, 1/3), col.to.transform = c(2, 5))
l.profile.mle <- profile.fwdtr(l.mle)
fwdtrProfilePlot(l.profile.mle)
}

\keyword{multivariate}
