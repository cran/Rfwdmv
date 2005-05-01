\name{threetwo.dat}
\alias{threetwo.dat}
\docType{data}

\title{Three clusters two outliers}

\description{
These data contain again the two clusters of the 60:80 data, but now with the addition of a third cluster, units 141-158 and two outliers, units 159 and 160.  The sizes of the groups are therefore 80, 60, 18 and 2
}

\usage{data(threetwo.dat)}
\format{
  A data frame with 160 observations on the following 2 variables.
  \describe{
    \item{y1}{a numeric vector.}
    \item{y2}{a numeric vector.}
  }
}

\details{
The scatter plot shows that the second compact cluster of 18 observations is near the longer axis of the dispersed cluster of 80.  The two outliers are together, approximately across the centroid of the dispersed group from the cluster of 60.
}

\source{
Atkinson, Riani and Cerioli (2004), p. 588-589; http://www.riani.it/arc.
}

\examples{
data(threetwo.dat)
}

\keyword{datasets}
