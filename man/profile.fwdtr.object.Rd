\name{profile.fwdtr.object}
\alias{profile.fwdtr.object}
\docType{data}
\title{profile.fwdtr.object}
\description{
An object containing the values of the profile log-likelihood for each transformation parameter
specified
in function \code{fwdtr}.  
The class attribute is set to \code{profile.fwdtr}.
}

\format{
 A \code{profile.fwdtr.object} is a list with the following elements:

\describe{
  
    \item{call}{the matched call.}
    \item{lambda}{a sequence of values from \code{bounds[1]} to \code{bounds[2]} with step 0.1.  
    This vector will contain the x-coordinates of the profile log-likelihoods}
    \item{profile}{a list containing the values of the profile log-likelihoods for each variable. 
    First element of the list is associated with first variable and so on.}
    \item{ci}{a list containing the lower and upper values of the confidence intervals
             for each variable. First element of the list is associated with first variable and so on.}
    \item{x.names}{a character vector containing the names of the columns which are investigated
          for transformation.}
    \item{step.fwd}{a scalar containing the step of the forward search for which 
                  the profile log-likelihoods are computed.}
             
    \item{onepar}{logical. If \code{TRUE} a common value \code{lambda} has been used for estimation}
    \item{p}{an integer, the number of variables in the data.}
    
    \item{Mle}{a vector containing the maximum likelihood estimates of the
               transformation parameters at step \code{step.fwd} of the forward search.}
 }

}
\details{
Class \code{profile.fwdtr} object is created by function \code{\link{profile.fwdtr}}
}
\seealso{

  \describe{
    \item{\code{\link{fwdtr}}}{computes maximum likelihood estimates of transformation parameters
    in each step of the forward search.} 
\item{\code{\link{fwdtrProfilePlot}}}{Plot profile loglikelihoods of transformation parameters in a selected step of the forward search }
} 
}

\examples{ 

data(mussels.dat)

## Forward search on untrasformed data 
## Compute max lik. estimates of tramsformation parameters

l.mle<-fwdtr(mussels.dat)

## Compute profile loglikelihoods for transformation parameters in the last step of the search 
## and create an object of class profile.fwdtr

l.profile.mle <- profile.fwdtr(l.mle)

## plot the profile loglikelihoods of transformation parameters for each variable

fwdtrProfilePlot(l.profile.mle)

}

\keyword{datasets}
