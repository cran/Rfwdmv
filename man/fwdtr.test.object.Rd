\name{fwdtr.test.object}
\alias{fwdtr.test.object}
\docType{data}
\title{fwdtr.test.object}
\description{
An object containing the values of the signed likelihood ratio test estimates of transformation parameters for selected variables and selected values of lambda as specified in function fwdtr.test.  The class attribute is set to \code{fwdtr.test}.
}

\format{

A \code{fwdtr.test.object} is a list with the following elements:

\describe{
  
    \item{call}{the matched call.}
    \item{test}{a list containing a collection of matrices. Each matrix has a number 
    of columns equal to the length of
    vector lambda.around. The first matrix contains the values
    of the signed sqrt likelihood ratio test for the first element in vector col.to.compare 
    supplied in function lambda.test.fwdmv.  More in detail, the first column contains the values of the signed 
    sqrt likelihood ratio test associated with first element of lambda.around. The second matrix 
    contains the values
    of the signed sqrt likelihood ratio test for the second element in vector col.to.compare, ecc... }
    \item{n}{an integer, the number of units in the data.}
    \item{m}{an integer, the number of units in the subset in the first step of the forward search.}
    \item{col.names}{a character vector containing the names of the columns of the data set.}
    \item{col.to.compare}{integer vector containing the numbers associated to 
    the variables for which the fan plot has been produced.}
    \item{lambda.around}{Vector specifying for which values of lambda the signed sqrt likelihood ratio has been computed.}
    \item{message}{Message that warn the user on the number of times that \code{optim} failed to converge properly. If \code{NULL} there are no convergence problem. See \link{optim} for details.}
   }
}


\details{
  Class \code{fwdtr.test} object is created by function \code{\link{fwdtr.test}}.
}

\seealso{
  \describe{
    \item{\code{\link{fwdtr}}}{computes maximum likelihood estimates of transformation parameters
    in each step of the forward search.} 
    }
}

\examples{

data(mussels.dat)

## reasonable values of transformation parameters for the data set
## found using procedure fwdtr

lambda.R <- c(0.5, 0, 0.5, 0, 0)
l.rat <- fwdtr.test(mussels.dat, lambda.R, col.to.compare = 1:5)

## l.rat is a fwdtr.test object

plot.fwdtr.test(l.rat)
}

\keyword{datasets}
