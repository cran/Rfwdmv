fwdmvCovariancePlot <- function(x, id = FALSE, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  p <- x$p
  initial <- x$initial
  Covariance <- x$Cov

  M <- m:n
	n.cols <- (p^2 + p) / 2

	if(is.character(id) && casefold(id) == "all") {
		id <- FALSE
		identify.all <- TRUE
	}
	
	else
	  identify.all <- FALSE

  n.groups <- length(Covariance)
  y.lim <- range(unlist(Covariance))

  ans <- list()

  if(id && !interactive())
	  return(ans)

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- "Elements of Covariance Matrix"
    x.label <- "Subset Size"
    main.label <- ""
  }

  if(id || identify.all)
    x.lim <- c(m, n + 0.1*(n-m))
  else
    x.lim <- c(m,n)

	var.cols <- 1:p
	var.cols <- (var.cols^2 + var.cols) / 2
	cov.cols <- setdiff(1:n.cols, var.cols)

  matplot(M,
					Covariance[[1]][, var.cols],
					type = "l",
					xlim = x.lim,
					ylim = y.lim,
					xlab = x.label,
					ylab = y.label,
					main = main.label,
					lty = 1,
					col = 1,
					lwd = 2)

  matlines(M,
					 Covariance[[1]][, cov.cols],
					 type = "l",
					 xlim = x.lim,
					 ylim = y.lim,
					 xlab = x.label,
					 ylab = y.label,
					 main = main.label,
					 lty = 1,
					 col = 1)

  if(n.groups > 1) {
    for(i in 2:n.groups) {
      matlines(M, Covariance[[i]][ , var.cols], lty = i, col = i, lwd = 2)
      matlines(M, Covariance[[i]][ , cov.cols], lty = i, col = i)
		}
  }
  
  if(!initial) {
    legend(m,
           y.lim[2],
           legend = x$group.names,
           col = 1:n.groups,
           lty = 1:n.groups)
  }

  if(id || identify.all) {
    xs <- matrix(rep(1:p, p), p, p)
    xs <- xs[row(xs) <= col(xs)]
    ys <- matrix(rep(1:p, each = p), p, p)
    ys <- ys[row(ys) <= col(ys)]
		diag.index <- xs == ys

    xs <- x$data.names[[2]][xs]
    ys <- x$data.names[[2]][ys]

		cov.labels <- character(n.cols)
		cov.labels[diag.index] <- paste(" var(", xs[diag.index], ")", sep = "")
    cov.labels[!diag.index] <- paste(" cov(", xs[!diag.index], ", ", ys[!diag.index], ")", sep = "")
	}

  if(id) {
    cat("\nLeft click to highlight trajectories. Right click\n")
    cat("to return to the R command prompt.\n\n")

    Cov <- Covariance[[1]]
    if(n.groups > 1) {
      for(i in 2:n.groups)
        Cov <- cbind(Cov, Covariance[[i]])
    }

    selected <- integer(0)
    while(length(selected) < dim(Cov)[2]) {
      xy <- locator(1)
      if(is.null(xy))
        break

      xy$x <- round(xy$x)
      if(xy$x > n)
        xy$x <- n
      xy$x <- xy$x - m + 1
      step.dist <- abs(Cov[xy$x, ] - xy$y)
      names(step.dist) <- 1:dim(Cov)[2]
      sel <- which(step.dist == min(step.dist))

      if(!is.element(sel, selected)) {
        selected <- c(selected, sel)
        color <- as.integer((sel - 1) / ((p^2 + p)/2)) + 1
        this <- sel - 1
        this <- this %% ((p^2 + p)/2) + 1
        text(n, Cov[n-m+1, sel], cov.labels[this], adj = 0, col = color)
      }
    }
    ans$selected <- selected
  }

  if(identify.all)
		for(i in 1:n.groups)
		  for(j in 1:n.cols)
		    text(n, Covariance[[i]][n-m+1, j], cov.labels[j], adj = 0, col = i)

  invisible(ans)
}


