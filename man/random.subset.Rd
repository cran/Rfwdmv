\name{random.subset}
\alias{random.subset}

\title{Random fwdmv Starting Subset}

\description{
Generates a random starting subset for multivariate forward searches.
}

\usage{
random.subset(X, n.bsb)
}

\arguments{
  \item{X}{a matrix or data frame containing the data.}
  \item{n.bsb}{a positive integer specifying the size of the subset.}
}

\value{
a numeric vector containing the random subset.
}

\seealso{\code{\link{fwdmv}}, \code{\link{fwdmv.init}}.}

\examples{
data(bank.dat)
random.subset(bank.dat, 7)

# Most common use: #
fwdmv(bank.dat, bsb = random.subset)
}

\keyword{multivariate}

