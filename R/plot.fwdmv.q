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


