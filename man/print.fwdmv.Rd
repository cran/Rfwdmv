\name{print.fwdmv}

\alias{print.fwdmv}

\title{Print method for fwdmv objects.}

\description{
Displays a brief summary of the forward search contained in an fwdmv object.
}

\usage{
\method{print}{fwdmv}(x, ...)
}

\arguments{
  \item{x}{an \code{fwdmv} object.}
  \item{\dots}{additional (and unused) arguments.}
}

\value{
\code{x} is invisibly returned.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv}}}

\examples{
  data(fondi.dat)
  fondi.1 <- fwdmv(fondi.dat)
  print(fondi.1)  ## Equivalent to fondi.1
}

\keyword{methods}
