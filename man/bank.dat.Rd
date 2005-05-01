\name{bank.dat}
\alias{bank.dat}
\docType{data}
\title{Swiss bank notes}

\description{
The data are readings on 200 Swiss bank notes, 100 of which are genuine and 100 forged.  However, the structure of the samples may not be quite that simple.  The forged notes have all been detected and withdrawn from circulation.  To provide a useful comparison, the genuine notes are likewise used.  So some of the notes in either group may have been misclassified.  A second complication is that there may be more than one forger at work.
}

\usage{data(bank.dat)}

\format{
  A data frame with 200 observations on the following 6 variables.
  \describe{
    \item{y1}{a numeric vector, the length of bank note near the top.}
    \item{y2}{a numeric vector, the left-hand height of bank note.}
    \item{y3}{a numeric vector, the right-hand height of bank note.}
    \item{y4}{a numeric vector, the distance from bottom of bank note to beginning of patterned border.}
    \item{y5}{a numeric vector, the distance from top of bank note to beginning of patterned border.}
    \item{y6}{a numeric vector, the diagonal distance.}
  }
}

\details{
The first three variables are measurements of paper size while the fourth and fifth variables are measurements from the edge of the paper to the printed area. Only \code{y6} is solely on the printed area.  It measures the diagonal distance across the frame of the central illustration.
}

\source{
Flury and Riedwyl (1988), pp. 4-8.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 562-566; http://www.riani.it/arc.
}

\examples{
data(bank.dat)
plot(bank.dat)
}

\keyword{datasets}

