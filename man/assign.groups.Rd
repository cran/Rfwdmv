\name{assign.groups}
\alias{assign.groups}

\title{Tentative Group Assignments}

\description{
This function is used to assign tentative groups to an \code{fwdmv} object containing and initial fit.
}

\usage{
assign.groups(object, groups)
}

\arguments{
  \item{object}{an \code{fwdmv} object containing and initial fit.}
  \item{groups}{a list containing the tentative groups.  Each group is given by a vector of positive integers and the groups must be mutually exclusive.}
}

\value{
an \code{fwdmv} object containing the tentative groups.
}

\author{Kjell Konis}

\seealso{\code{\link{fwdmv}}}

\examples{
data(fondi.dat)

g1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 53, 55, 56)

g2 <- c(57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)

# Create an fwdmv object containing just an initial fit #

fondi.init <- fwdmv(fondi.dat)

# Assign tentative groups #

fondi.tgs <- assign.groups(fondi.init, groups = list(g1, g2))
}

\keyword{utilities}
