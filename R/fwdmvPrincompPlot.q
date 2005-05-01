fwdmvPrincompPlot <- function(x, psfrag.labels = FALSE)
{  
  m <- x$m
  n <- x$n
  M <- m:n
  group.names <- x$group.names

  if(length(group.names) == 1 && group.names == "Unassigned")
    group.names <- ""

  ans <- list()

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- "Principal Components"
    x.label <- "Subset Size"
    main.label <- ""
  }

  eigenvalues <- eigenvalues.fwdmv(x)
  n.groups <- length(eigenvalues)
	for(i in 1:n.groups)
	  eigenvalues[[i]] <- eigenvalues[[i]] / apply(eigenvalues[[i]], 1, sum)

  matplot(M,
          eigenvalues[[1]],
          type = "l",
          ylim = c(0.0, 1.0),
          xlab = x.label,
          ylab = y.label,
          main = main.label,
          lty = 1,
          col = 1)

  if(n.groups > 1)
    for(i in 2:n.groups)
      matlines(M, eigenvalues[[i]], lty = i, col = i)

  if(nchar(group.names[1]))
    legend(m,
           1.0,
           legend = group.names,
           col = 1:n.groups,
           lty = 1:n.groups)

  invisible(ans)
}


