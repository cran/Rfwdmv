\name{fondi.dat}
\alias{fondi.dat}
\docType{data}

\title{Investment Funds Data}

\description{
These data contain information on 103 investment funds operating in Italy since April 1996.
}

\usage{data(fondi.dat)}

\format{
  A data frame with 103 observations on the following 3 variables.
  \describe{
    \item{y1}{a numeric vector, the short term (12 months) performance.}
    \item{y2}{a numeric vector, the medium term (36 months) performance.}
    \item{y3}{a numeric vector, the medium term (36 months) volatility.}
  }
}

\details{
The scatter plot matrix of the data shows there seem to be two clusters, with a few observations in between.
}

\source{
Zani (2000), p. 194.  An introduction to the data in English is given by Cerioli and Zani (2001).
}

\references{
Atkinson, Riani and Cerioli (2004), p. 592-593; http://www.riani.it/arc.
}

\examples{
data(fondi.dat)
}

\keyword{datasets}
