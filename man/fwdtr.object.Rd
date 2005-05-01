\name{fwdtr.object}
\alias{fwdtr.object}
\docType{data}
\title{fwdtr.object}
\description{
An object containing estimates of transformation parameters on multivariate data.  The class attribute is set to \code{fwdtr}.
}

\format{
A \code{fwdtr.object} is a list with the following elements:

\describe{
    \item{call}{the matched call.}

    \item{Unit}{a list containing the units forming the subset in each step of the search.}

    \item{n}{an integer, the number of units in the data.}

    \item{p}{an integer, the number of variables in the data.}

    \item{m}{an integer, the number of units in the initial subset.}

    \item{data}{a numeric matrix containing the data, the dimnames attribute is set to \code{NULL}.}

    \item{data.H0}{a numeric matrix containing the data transformed using lambda under the null hypotehsis.}

    \item{data.name}{the name of the data frame or matrix containing the multivariate data set.}

    \item{data.names}{a list of character vectors containing the row and column names of the data.}

    \item{Mle}{a matrix containing the estimates of transformation parameters in each step of the forward search.  The number of columns of matrix Mle is equal to the length of vector \code{col.to.transform} or is equal to 1 if \code{one.lambda = TRUE}.}

    \item{H0}{k x 4 matrix. The first column contains the values of lambda (transformation paramters) which have been used to transform the variables to order Mahalanobis distances during the forward search. The second column contains the integers associated to the columns for which we computed maximum likelihood estimates. The third and fourth column contain respectively the constrained lower and upper bounds of the estimates of the values of the transformation parameters.}

    \item{forced.onepar}{Logical value. if \code{TRUE} a common value of lambda is estimated for all variables specified in \code{col.to.transform}.}

    \item{dof}{scalar integer value containing the degrees of freedom to be used in the likelihood ratio test.}
    }

}
\details{
Class \code{fwdtr} object is created by function \code{\link{fwdtr}}.
}

\seealso{
  \describe{
    \item{\code{\link{fwdmv.init}}}{fit an ititial multivariate forward search.}

    \item{\code{\link{fwdmv}}}{fit a multivariate forward search with user specified tentative groups.}
  }
}

\examples{
data(mussels.dat)
l.mle <- fwdtr(mussels.dat)

#l.mle is a fwdtr object.
}

\keyword{datasets}
