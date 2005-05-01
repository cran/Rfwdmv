fwdmvEigenvectorPlot <- function(x, which.vector = 1, correlation = FALSE,
                                   psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  p <- x$p
  M <- m:n

  ans <- list()

  group.names <- x$group.names
  data.names <- x$data.names

  if(length(group.names) == 1 && group.names[1] == "Unassigned")
    group.names <- ""

  if(p <= 3)
    lty <- c(1,2,4)[1:p]
  else
    lty <- 1:p

  eigenvectors <- eigenvectors.fwdmv(x, which.vector = which.vector)
  n.groups <- length(eigenvectors)

	if(correlation) {
	  eigenvalues <- eigenvalues.fwdmv(x)
    for(i in 1:n.groups)
		  eigenvectors[[i]] <- eigenvectors[[i]] / sqrt(eigenvalues[[i]][, which.vector])
	}

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    x.label <- "Subset Size"
    y.label <- paste("Elements of Eigenvector", which.vector)
    main.label <- ""
  }

  old.par <- par(pty = "s")
  on.exit(par(old.par))

  matplot(M,
          eigenvectors[[1]],
          type = "l",
          ylim = c(-1.0, 1.0),
          xlab = x.label,
          ylab = y.label,
          main = main.label,
          col = 1,
          lty = lty)

  if(n.groups > 1)
    for(i in 2:n.groups)
      matlines(M,
               eigenvectors[[i]],
               col = i,
               lty = lty)

  curve.names <- apply(expand.grid(group.names, data.names[[2]]), 1,
                       paste, collapse = ": ")

  if(min(eigenvectors[[1]][1,]) + 1.0 > 1.0 - max(eigenvectors[[1]][1,])) {
    y.pos <- -1.0
    yjust <- 0
  }
  else {
    y.pos <- 1.0
    yjust <- 1
  }

  legend(m,
         y.pos,
         legend = curve.names,
         col = rep(1:n.groups, each = p),
         lty = rep(lty, n.groups),
         yjust = yjust,
				 ncol = n.groups)

  invisible(ans)
}


