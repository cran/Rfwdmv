\name{FScomputeEnvelope}
\alias{FScomputeEnvelope}

\title{Compute a Confidence Envelope}

\description{
Computes a confidence envelope for the Random Starts forward plot.
}

\usage{
FScomputeEnvelope(n, p, m, probs = c(0.01, 0.5, 0.99), tolerance = 0.001)
}

\arguments{
\item{n}{an integer value giving the number of observations in the data.}

\item{p}{an integer value giving the dimension of the data.}

\item{m}{an integer giving the size of the smallest subset in the forward search.}

\item{probs}{a numeric vector containing the probabilities used to compute the quantiles.}

\item{tolerance}{a positive numeric value giving the convergence tolerance.}
}

\value{
a numeric matrix whose columns contain the quantiles of the given order statisics.  The row names contain the subset sizes and the columns names contain the quantiles specified in \code{probs}.
}

\seealso{
\code{\link{fwdmvRandomStart}}, 
\code{\link{plot.fwdmvRandomStart}}.
}

\examples{
X <- matrix(rnorm(1200), nrow = 200, ncol = 6)
temp <- fwdmvRandomStart(X)

# the envelope is computed by FScomputeEnvelope #
plot(temp, envelope = FALSE)
envelope <- FScomputeEnvelope(200, 6, 12)
domain <- as.numeric(row.names(envelope))
correction <- 200*6*(domain-2)/((200-1)*(domain-6-1))
envelope <- sqrt(sweep(envelope, 1, correction, "*"))
matlines(domain, envelope, lty = 2, col = 1)
}

\keyword{multivariate}

