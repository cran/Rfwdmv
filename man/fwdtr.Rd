\name{fwdtr}
\alias{fwdtr}

\title{ Maximum likelihood estimates of transformation parameters }
\description{This function computes maximum likelihood estimates of transformation parameters. It uses the multivariate version of the parametric family of power transformations introduced by Box and Cox (1964).
}
\usage{
fwdtr(X, bsb = ellipse.subset, n.bsb, lambda = 1, one.lambda = FALSE, col.to.transform = "all", boundaries = c(-3, 3))
}

\arguments{
  \item{X}{ a matrix or data frame containing a multivariate data set.}
  \item{bsb}{usually a function of two variables: a matrix \code{X} containing the multivariate data and the number of units in the initial subset \code{n.bsb}.  Alternatively, the initial subset my be specified directly by providing an integer vector containing the indices of the units to include in the initial subset.}
  \item{n.bsb}{the percentage of units forming the initial subset.  For example \code{n.bsb = 40} implies that we start the search using \code{as.integer(nrow(X) * 0.4)} units. The default is \code{n.bsb = 50}.}
  \item{lambda}{a scalar or a \code{k x 1} vector containing set of transformation parameters.  The ordering of Mahalanobis distances at each step of the forward search uses variables transformed with \code{lambda}. If \code{lambda} is a scalar all the variables in \code{col.to.transform} are transformed, for ordering Mahalanobis distances, using the commom supplied value. }
  \item{one.lambda}{ a logical value.  If \code{TRUE} a common value \code{lambda} is estimated for all variables specified in \code{col.to.transform}.}
  \item{col.to.transform}{a \code{k x 1} integer vector specifying the variables which must be transformed. If \code{col.to.transform = "all"} all variables (columns of matrix \code{X}) are considered for transformation.}
  \item{boundaries}{the upper and lower bounds for the estimates of the values of the 
  transformation parameters.}
}
\details{
 The analysis of data can often be improved by using transformed variables rather than the original variables themselves.  There are physical reasons why a transformation might be expected to be helpful in some examples. If the data arise from a counting process, they often have a Poisson distribution and the square root transformation will provide observations with an approximately constant variance, independent of the mean.  Similarly, concentrations are nonnegative variables and so cannot strictly be subject to additive errors of constant variance.  Unfortunately the estimated transformation and related test statistics  may be sensitive to the presence of one, or several, outliers.  With this function we use the forward search to see how  estimates of the transformation parameters evolve as we move through the ordered data. If a correct value of \code{lambda} has been found the parameter estimates will be stable until near the end of the search, where any outliers start to enter.
}
\value{
 a list with class \code{fwdtr}.
}
\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York. }
\author{ Fabrizio Laurini}


\seealso{\code{\link{fwdtr.object}}, \code{\link{profile.fwdtr.object}} }
 
\examples{ 
## Forward search on untrasformed data 
data(mussels.dat)
l.mle <- fwdtr(mussels.dat)
fwdtrMlePlot(l.mle)

## Forward search on transformed data as specified in vector lamdba 
data(mussels.dat)
l.mle <- fwdtr(mussels.dat,lambda = c(1, 0.5, 1, 0, 1/3))
fwdtrMlePlot(l.mle)

## estimate a common value of lambda for all the variables and use 1/3 to order Mahalanobis distances in each step of the search
data(mussels.dat)
l.mle <- fwdtr(mussels.dat, lambda = 1/3, one.lambda = TRUE)
fwdtrMlePlot(l.mle)

### Test variables 2 and 5
### The forward is based on untransformed data for variables 1, 3 and 4
### sqrt for variable 2 and third root for variable 5
l.mle <- fwdtr(mussels.dat, lambda = c(0.5, 1/3), col.to.transform = c(2, 5))
fwdtrMlePlot(l.mle)
}

\keyword{multivariate}
