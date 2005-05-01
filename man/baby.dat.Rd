\name{baby.dat}
\alias{baby.dat}
\docType{data}

\title{Babyfood data}

\description{
The responses are the density in centipoises at the time of manufacture (y1) and when measured 3 (y2), 6 (y3) and 9 months later (y4).  There are five explanatory variables.
}

\usage{data(baby.dat)}

\format{
  A data frame with 27 observations on the following 9 variables.
  \describe{
    \item{x1}{a numeric vector.}
    \item{x2}{a numeric vector.}
    \item{x3}{a numeric vector.}
    \item{x4}{a numeric vector.}
    \item{x5}{a numeric vector.}
    \item{y1}{a numeric vector, the initial viscosity of the babyfood.}
    \item{y2}{a numeric vector, the viscosity of the babyfood after three months storage.}
    \item{y3}{a numeric vector, the viscosity of the babyfood after six months storage.}
    \item{y4}{a numeric vector, the viscosity of the babyfood after nine months storage.}
  }
}

\details{
Box and draper (1987, p. 572) find  a linear model with terms x2, x3 and x5, as well as surprisingly the interaction x3:x4 in the absence of x4.  This model was suggested for all four responses.  It is generally agreed that such models, violating a marginality constraint, are undesirable: if the variables in this model are rescaled, the model will apparently change, a term in x3 appearing.
}

\source{
Box and Draper (1987), p. 265 present part of a larger data set on the storage of a babyfood. 
}

\references{
Atkinson, Riani and Cerioli (2004), p. 567; http://www.riani.it/arc.
}

\keyword{datasets}

