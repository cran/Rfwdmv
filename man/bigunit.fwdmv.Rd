\name{bigunit.fwdmv}
\alias{bigunit.fwdmv}

\title{Generate the Big Unit Matrix}

\description{
Returns a logical matrix with one column for each step in the Forward Search and one row for each unit in the data.  The (i,j) element is \code{TRUE} if unit i is in the subset during step j and \code{FALSE} otherwise.
}

\usage{
bigunit.fwdmv(x)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
}

\value{
a logical matrix.  See the description.
}

\references{Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. Springer-Verlag New York.}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv.object}}, \code{\link{fwdmv}}}

\examples{
data(fondi.dat)
fondi.1 <- fwdmv(fondi.dat)
fondi.bigUnit <- bigunit.fwdmv(fondi.1)
}

\keyword{utilities}

