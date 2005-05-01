\name{fwdmv.init}
\alias{fwdmv.init}

\title{Multivariate Forward Search for Ungrouped Data}

\description{
This function computes a multivariate forward search for ungrouped data.  Several diagnostic statistics are monitored during the search: see \code{fwdmv.object}.  Note that this function is called by \code{\link{fwdmv}} when no tentative groups are specified.  It is recommended that \code{\link{fwdmv}} be used for all multivariate forward searches.
}

\usage{
fwdmv.init(X, bsb = ellipse.subset, scaled = TRUE)
}

\arguments{
  \item{X}{a matrix or data frame containing the multivariate data set.}
  \item{bsb}{a function of two variables: the multivariate data in matrix form \code{X} 
    and the number of units in the initial subset \code{size}. 
    If bsb= bb.subset the initial susbet is found using robust bivariate boxplots. The default is to 
    use robustly centered ellipses  \code{ellipse.subset}.
    Alternatively, the initial subset my be specified directly by providing an integer vector 
    containing the indices of the units to be in the initial subset.}
  \item{scaled}{a logical value.  If \code{TRUE} then scaled Mahalanobis distances are used during the forward search.}
}

\details{
This function computes the Forward Search described in chapter 3 of ARC.  The initial subset can be 
specified directly in the argument \code{bsb} or computed from the data.  
By default \code{bsb} is a function for computing the initial subset using robustly centered ellipses.  
Given a subset of \code{m} units the next subset is the \code{m+1} 
units with smallest Mahalanobis distances calculated using the center and covariance matrix 
estimates of the units currently in the subset.  
This process is repeated until the subset contains all of the units and several diagnostic 
statistics are computed for each subset.
}

\value{
  aa \code{fwdmv} object.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}}

\examples{
data(fondi.dat)
fondi.init <- fwdmv.init(fondi.dat)

data(fondi.dat)
#### find the intial subset using robust bivariate ellipses  
#### start with an initial subset size of 17 units
fondi.init <- fwdmv.init(fondi.dat,bsb=ellipse.subset(fondi.dat,17))

data(fondi.dat)
#### find the intial subset using robust bivariate boxplots and 
#### start with an initial subset size of 17 units and 
fondi.init <- fwdmv.init(fondi.dat,bsb=bb.subset(fondi.dat,17))

}

\keyword{multivariate}
