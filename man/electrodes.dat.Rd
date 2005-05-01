\name{electrodes.dat}
\alias{electrodes.dat}
\docType{data}

\title{Electrodes data}

\description{
The data are measurements from machines manufacturing supposedly identical electrodes.  The electrodes are shaped rather like nipples.  There are five measurements on each: \code{y1}, \code{y2} and \code{y5} are diameters, while \code{y3} and \code{y4} are lengths.  Fifty electrodes from each machine have been measured.
}

\usage{data(electrodes.dat)}

\format{
  A data frame with 100 observations on the following 5 variables.
  \describe{
    \item{y1}{a numeric vector, the first diameter.}
    \item{y2}{a numeric vector, the second diameter.}
    \item{y3}{a numeric vector, the first length.}
    \item{y4}{a numeric vector, the second length.}
    \item{y5}{a numeric vector, the third diameter.}
  }
}

\details{
For reasons of commercial secrecy, the data have been transformed by subtracting constants from the variables.
}

\source{
The data come from an unpublished University Ph.D. thesis by Kreuter.  The data are given and described by Flury and Riedwyl (1998), pp. 128-132.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 579-580; http://www.riani.it/arc.
}

\examples{
data(electrodes.dat)
}

\keyword{datasets}
