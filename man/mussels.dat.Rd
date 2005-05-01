\name{mussels.dat}
\alias{mussels.dat}
\docType{data}

\title{Horse mussels}

\description{
These data contain 82 observations on Horse mussels from New Zealand.
}

\usage{data(mussels.dat)}

\format{
  A data frame with 82 observations on the following 5 variables.
  \describe{
    \item{width}{a numeric vector, the shell length in mm.}
    \item{height}{a numeric vector, the shell width in mm.}
    \item{length}{a numeric vector the shell height in mm.}
    \item{shell}{a numeric vector the shell mass in grams.}
    \item{mass}{a numeric vector, the muscle mass in grams.}
  }
}

\details{
In this example the effect of outliers is masked unless the search uses a suitable transformation.
}

\source{
The data were introduced by Cook and Weisberg (1994), p. 161 who treat them as regression with muscle mass, the edible portion of the mussel, as response.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 568-569; http://www.riani.it/arc.
}

\examples{
data(mussels.dat)
}

\keyword{datasets}
