\name{milk.dat}
\alias{milk.dat}
\docType{data}

\title{Milk data}

\description{
These data contain measurements on the composition of 85 samples of milk.
}

\usage{data(milk.dat)}

\format{
  A data frame with 85 observations on the following 8 variables.
  \describe{
    \item{y1}{a numeric vector, the density.}
    \item{y2}{a numeric vector, the fat content in grams/litre.}
    \item{y3}{a numeric vector, the protein content in grams/litre.}
    \item{y4}{a numeric vector, the casein content in grams/litre.}
    \item{y5}{a numeric vector, the cheese dry substance measured in the factory in grams/litre.}
    \item{y6}{a numeric vector, the cheese dry substance measured in the laboratory in grams/litre.}
    \item{y7}{a numeric vector, the milk dry substance in grams/litre.}
    \item{y8}{a numeric vector, the cheese produced in grams/litre.}
  }
}

\details{
The scatter plot matrix shows in several panels a strong rising diagonal structure.
}

\source{
Daudin, Duby and Trecourt (1988).
}

\references{
Atkinson, Riani and Cerioli (2004), p. 571-572; http://www.riani.it/arc.
}

\examples{
data(milk.dat)
pairs(milk.dat)
}

\keyword{datasets}
