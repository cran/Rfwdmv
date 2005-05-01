\name{fwdtr.test}
\alias{fwdtr.test}

\title{Multivariate Fan Plot}
\description{
Confirmatory signed square root likelihood ratio tests of a suggested transformation for specified variables around a series of values of lambda.
}
\usage{
fwdtr.test(X, parameters,  n.bsb = 50, col.to.compare = "all", lambda.around = c(-1, -0.5, 0, 0.5, 1), one.lambda = FALSE)
}

\arguments{
  \item{X}{a matrix of dimension n x p, or data frame containing a multivariate data set. }
  \item{parameters}{a vector of length p=ncol(X) specifying a reasonable set of transformations for the columns of the multivariate data set.}
  \item{n.bsb}{the percentage of units forming the initial subset.  For example \code{n.bsb = 40} implies that we start the search using \code{as.integer(nrow(X)*0.4)} units.  The default is \code{n.bsb=50}.}
  \item{col.to.compare}{a \code{k x 1} integer vector specifying the variables for which likelihood ratio tests have to be produced.  For example, if \code{col.to.compare = c(2, 4)}, the signed likelihood ratio tests are produced for the second and the fourth column of matrix X.  If \code{col.to.compare = "all"} the all variables (columns of matrix \code{X}) are considered.}
  \item{lambda.around}{a numeric vector specifying for which values of lambda to compute the likelihood ratio test.  If this argument is omitted, the function produces for each variable specified in \code{col.to.compare}  the likelihood ratio tests associated to the five most common values of lambda \code{(-1, -0.5, 0, 0.5, 1)}.}
  \item{one.lambda}{a logical value.  If \code{TRUE} a common value \code{lambda} is tested for all variables specified in \code{col.to.transform}.}
}
\details{This function produces confirmatory tests of a suggested transformation.  We expand each transformation parameter
in turn around the five common values of lambda (-1, -0.5, 0, 0.5, 1), using the values of the vector \code{parameters} for transforming the remaining variables of the data set.  In this way we turn a multivariate problem into a series of univariate ones.  In each search we can test the transformation by comparing the likelihood ratio test with a chisquare on 1 degree of freedom.  We use the signed square root of the likelihood ratio test in order to learn whether lower or higher values of lambda are indicated. The plot is thus a version of the fan plot for multivariate data.
}
\value{ a list with class \code{fwdtr.test}.
}
\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004)  Exploring Multivariate Data with the Forward Search.  Springer-Verlag New York.}
\author{Fabrizio Laurini}


\seealso{\code{\link{fwdtr.object}}, \code{\link{profile.fwdtr.object}}}
\examples{
data(mussels.dat)

## reasonable values of transformation parameters for the data set

lambda.R <- c(0.5, 0, 0.5, 0, 0)
l.rat <- fwdtr.test(mussels.dat, lambda.R)
plot.fwdtr.test(l.rat)

## Produce a fan plot for columns 2 and 4 of dataset mussels.dat

l.rat <- fwdtr.test(mussels.dat, parameters = lambda.R, col.to.compare = c(2,4))
plot.fwdtr.test(l.rat)

## reasonable values of transformation parameters for the data set

lambda.R <- c(0.5,0,0.5,0,0)
lambda.around <- c(0,1/3,0.5)

## Produce a fan plot for column 2 of dataset mussels.dat
## The values of lambda which are tested are log, third root
## and square root

l.rat <- fwdtr.test(mussels.dat, parameters = lambda.R, col.to.compare = 2, lambda.around = lambda.around)
plot.fwdtr.test(l.rat)
}

\keyword{multivariate}
