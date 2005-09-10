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


