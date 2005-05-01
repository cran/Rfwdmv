\name{sixtyeighty.dat}
\alias{sixtyeighty.dat}
\docType{data}

\title{60:80 data}

\description{
An example with two clusters.  The first 80 units form a rather diffuse group, while the remaining 60 units form a tight cluster.
}

\usage{data(sixtyeighty.dat)}

\format{
  A data frame with 140 observations on the following 2 variables.
  \describe{
    \item{y1}{a numeric vector.}
    \item{y2}{a numeric vector.}
  }
}

\details{
When we fit a single model to these data, neither standard classical methods nor very robust methods yield Mahalanobis distances which unambiguously show that there are two different groups of observations.
}

\source{
Atkinson, Riani and Cerioli (2004), p. 586-587; http://www.riani.it/arc.
}

\examples{
data(sixtyeighty.dat)
}

\keyword{datasets}
