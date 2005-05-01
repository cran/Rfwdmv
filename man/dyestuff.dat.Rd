\name{dyestuff.dat}
\alias{dyestuff.dat}
\docType{data}

\title{Dyestuff data}

\description{
The data arise in a study of the manufacture of a dyestuff.  There are 64 observations at the points of a \code{2^6} factorial and three responses: strength, hue and brightness.
}

\usage{data(dyestuff.dat)}

\format{
  A data frame with 64 observations on the following 6 variables.
  \describe{
    \item{x1}{a numeric vector.}
    \item{x2}{a numeric vector.}
    \item{x3}{a numeric vector.}
    \item{y1}{a numeric vector giving the strength.}
    \item{y2}{a numeric vector giving the hue.}
    \item{y3}{a numeric vector giving the brightness.}
  }
}

\details{
The original data set contained 6 variables.  However only 3 of the six explanatory variables have a significant effect on the three responses.  Here we report only the three explanatory variables which have been found to be significant.
}

\source{
Box and Draper (1987), pp. 114-115.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 570; http://www.riani.it/arc.
}

\examples{
data(dyestuff.dat)
}

\keyword{datasets}
