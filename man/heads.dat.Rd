\name{heads.dat}
\alias{heads.dat}
\docType{data}

\title{Swiss heads data}

\description{
These data contain six dimensions in millimetres of the heads of 200 Swiss soldiers.
}

\usage{data(heads.dat)}

\format{
  A data frame with 200 observations on the following 6 variables.
  \describe{
    \item{y1}{a numeric vector, the minimal frontal breadth.}
    \item{y2}{a numeric vector, the breadth of angulus mandibulae.}
    \item{y3}{a numeric vector, the true facial height.}
    \item{y4}{a numeric vector, the length from glabella to apex nasi.}
    \item{y5}{a numeric vector, the length from tragion to nasion.}
    \item{y6}{a numeric vector, the length from tragion to gnathion.}
  }
}

\details{
The data were collected to determine the variability in size and shape of heads of young men in order to help in the design of a new protection mask for the Swiss army.
}

\source{
The data are described by Flury and Riedwyl (1998), p. 218 and also by Flury (1997), p. 6.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 592-593; http://www.riani.it/arc.
}

\examples{
data(heads.dat)
}

\keyword{datasets}
