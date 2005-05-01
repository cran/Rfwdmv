\name{quality.dat}
\alias{quality.dat}
\docType{data}

\title{Quality of life data}

\description{
Indices of the quality of life in the provinces of Italy. 
}

\usage{data(quality.dat)}

\format{
  A data frame with 103 observations on the following 6 variables.
  \describe{
    \item{y1}{a numeric vector, the average amount of bank deposits per inhabitant.}
    \item{y2}{a numeric vector, the number of robberies per 100,000 inhabitants.}
    \item{y3}{a numeric vector, the number of housebreakings per 100,000 inhabitants.}
    \item{y4}{a numeric vector, the number of suicides, committed or attempted, per 100,000 inhabitants.}
    \item{y5}{a numeric vector, the number of gyms per 100,000 inhabitants.}
    \item{y6}{a numeric vector, the average expenditure on theatre and concerts per inhabitant.}
  }
}

\details{
The data are not provided on the original variables listed above, but rather on a scaled version of them.  Scaling is performed by dividing each response value by the maximum reading for that response.
}

\source{
These data have been published on the Italian financial newspaper Il Sole 24 ore 2001.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 573-575; http://www.riani.it/arc.
}

\examples{
data(quality.dat)
}

\keyword{datasets}
