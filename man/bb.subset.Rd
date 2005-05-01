\name{bb.subset}
\alias{bb.subset}

\title{Initial Subset by Bivariate Box Plots}

\description{
Computes the initial subset using robust bivariate boxplots as described in ARC. 
}

\usage{
bb.subset(X, size)
}

\arguments{
  \item{X}{a matrix or data frame containing a multivariate data set.} 
  \item{size}{the size of the initial subset.}
}

\value{
An integer vector containing the indexes of the units in the initial subset.
}

\references{
Zani, S., Riani M. and Corbellini A. (1998), 
"Robust Bivariate Boxplots and Multiple Outlier Detection", 
Computational Statistics and Data Analysis (1998), pp. 257-270.

Atkinson, A. C., Riani, M. and Cerioli, A. (2004) Exploring Multivariate Data with the Forward Search. 
Springer-Verlag New York.}

\author{Aldo Corbellini}

\seealso{ \code{\link{fwdmv}}   \code{\link{fwdmv.init}}  \code{\link{fwdmvPrePlot}}  }

\examples{
data(fondi.dat)
fondi.fwdmv <- fwdmv(fondi.dat, bsb = bb.subset)

### start with a subset size m=22
data(fondi.dat)
fondi.fwdmv <- fwdmv(fondi.dat, bsb = bb.subset(fondi.dat,22))
}

\keyword{utilities}

