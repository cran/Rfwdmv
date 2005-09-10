\name{plot.fwdmvRandomStart}
\alias{plot.fwdmvRandomStart}

\title{Redraw the Random Start Forward Plot}

\description{
Redraws the random start forward plot.
}

\usage{
\method{plot}{fwdmvRandomStart}(x, envelope = TRUE, ...)
}

\arguments{
  \item{x}{an \code{fwdmvRandomStart} object created by \code{\link{fwdmvRandomStart}}.}
  \item{envelope}{a logical flag.  If \code{TRUE} then a confidence envelope is included in the plot.}
  \item{\dots}{additional argument required by the geneic plot function.}
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
  \code{x} is invisibly returned.
}

\seealso{\code{\link{fwdmvRandomStart}}}

\examples{
data(bank.dat)
bank.rfs <- fwdmvRandomStart(bank.dat, plot.it = FALSE)
plot(bank.rfs)
# use identify on bank.rfs #
}

\keyword{hplot}

