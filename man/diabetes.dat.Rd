\name{diabetes.dat}
\alias{diabetes.dat}
\docType{data} 

\title{Diabetes data}

\description{
This data set consists of 145 observations on diabetes patients.  These data have been used in the statistical literature as a difficult example of cluster analysis.
}

\usage{data(diabetes.dat)}

\format{
  A data frame with 145 observations on the following 3 variables.
  \describe{
    \item{y1}{a numeric vector, the plasma glucose response to oral glucose.}
    \item{y2}{a numeric vector, the plasma insulin response to oral glucose.}
    \item{y3}{a numeric vector, the degree of insulin resistance.}
  }
}

\details{ 
\code{y1} and \code{y2} are responses to oral glucose, y3 is insulin resistance.  The scatter plot matrix of the data shows that there seems to be a central cluster and two "arms" forming separate clusters.
}

\source{
Reaven and Miller (1979), Fraley and Raftery (1998).
}

\references{
Atkinson, Riani and Cerioli (2004), p. 594-595; http://www.riani.it/arc.
}

\keyword{datasets}
