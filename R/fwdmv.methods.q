# S style OOP methods #

print.fwdmv <- function(x, ...)
{
  groups <- x$groups
  m <- x$m
  n <- x$n
  data.name <- x$data.name
  intial <- x$initial
  
  if(n.groups <- length(groups)) {
    cat(paste("\n\tA multivariate forward search on the data ", data.name, ".\n", sep = ""))
    cat(paste("\tThere are ", n.groups, " tentative group(s).\n", sep = ""))
    cat(paste("\tThe initial subset contains ", m, " units and the final subset ", n, " units.\n\n", sep = ""))
  }

  else {
    cat(paste("\n\tA multivariate forward search on the data ", data.name, ".\n", sep = ""))
    cat("\tA preliminary forward search.\n")
    cat(paste("\tThe initial subset contains ", m, " units and the final subset ", n, " units.\n\n", sep = ""))
  }

  invisible(x)
}


plot.fwdmv <- function(x, ...)
{
	if(!interactive())
	  return()

  choices <- c(" Change Plot",
               " Covariance Plot",
               " Determinant Plot",
               " Distance Plot",
               " Principal Components Plot",
               " Entry Order Plot",
               " Gap Plot",
               " Minimum and Maximum Distances Plot",
               " Pairs Plot",
							 " Eccentricity Plot",
							 " Eigenvector Plot")

  choice <- menu(choices, title = "Select a Plot")

  ans <- switch(choice,
    fwdmvChangePlot(x, ...),
    fwdmvCovariancePlot(x, ...),
    fwdmvDeterminantPlot(x, ...),
    fwdmvDistancePlot(x, ...),
    fwdmvPrincompPlot(x, ...),
    fwdmvEntryPlot(x, ...),
    fwdmvGapPlot(x, ...),
    fwdmvMinmaxPlot(x, ...),
    fwdmvPairsPlot(x, ...),
    fwdmvEccentricityPlot(x, ...),
    fwdmvEigenvectorPlot(x, ...))

  invisible(ans)
}


