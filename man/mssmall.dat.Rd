\name{mssmall.dat}
\alias{mssmall.dat}
\docType{data}

\title{Muscular dystrophy data (small)}

\description{
These data  are referred to Duchenne muscular dystrophy (DMD), a genetically transmitted disease passed from a mother to her children.  Affected male offspring may unknowingly carry the disease but female offspring with the disease die at a young age.  Although carriers of DMD usually have no physical symptoms, they tend to exhibit elevated levels of serum markers.  In addition the levels of these enzymes may also depend on age and season.  Levels of the enzymes were measured in non carriers and in a group of carriers using standard laboratory procedures.
}

\usage{data(mssmall.dat)}

\format{
  A data frame with 73 observations on the following 6 variables.
  \describe{
    \item{y1}{a numeric vector, the age.}
    \item{y2}{a numeric vector, the month of the year.}
    \item{y3}{a numeric vector, the level of creatine kinase.}
    \item{y4}{a numeric vector, the level of hemopexin.}
    \item{y5}{a numeric vector, the level of lactate dehydrogenase.}
    \item{y6}{a numeric vector, the level of pyruvate kinase.}
  }
}

\details{
The first two serum markers, \code{y3} and \code{y4}, may be measured rather inexpensively from frozen serum.  The second two, \code{y5} and \code{y6}, require fresh serum.  An important scientific problem is whether use of the expensive second pair of readings causes an appreciable increase in the detection rate.
}

\source{
Rencher (1995), p. 170.
}

\references{
Atkinson, Riani and Cerioli (2004), p. 581-585; http://www.riani.it/arc}
}

\examples{
data(mssmall.dat)
}

\keyword{datasets}
