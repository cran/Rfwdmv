fwdmvDeterminantPlot <- function(x, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  M <- m:n
  group.names <- x$group.names

  ans <- list()

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- "Determinant"
    x.label <- "Subset Size"
    main.label <- ""
  }

  determinants <- x$Determinant
  n.groups <- length(determinants)
  
  y.lim <- c(0, max(unlist(determinants)))

  if(n.groups > 1) {
    plot(M,
         determinants[[1]],
         type = "l",
         ylim = y.lim,
         xlab = x.label,
         ylab = y.label,
         main = main.label)

    for(i in 2:n.groups)
      lines(M, determinants[[i]], lty = i, col = i)
  }
  
  else
    plot(M,
         determinants[[1]],
         type = "l",
         xlab = x.label,
         ylab = y.label,
         main = main.label)

  if(group.names[1] != "Unassigned")
    legend(m,
           y.lim[2],
           legend = group.names,
           col = 1:n.groups,
           lty = 1:n.groups)

  invisible(ans)
}


