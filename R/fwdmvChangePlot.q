fwdmvChangePlot <- function(x, psfrag.labels = FALSE)
{
  n <- x$n
  m <- x$m
	groups <- x$groups
	group.names <- x$group.names
	unassigned <- x$unassigned
  distances <- x$Distances

	n.groups <- length(groups)
	n.unassigned <- length(unassigned)
  M <- (m+1):n
  ans <- list()
	plot.colors <- colors()[c(1, 260, 387, 490, 536, 640, 498, 552, 652)]
  
  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- "Positive Change in Distance"
    x.label <- "Subset Size"
    main.label <- ""
  }

  zz <- distances[, -1] > distances[, -(n-m+1)]
	
	if(length(groups)) {
	
		zz <- zz[c(unlist(groups), unassigned), ]
		zz <- zz * rep(1:(n.groups + 1), times = c(sapply(groups, length), n.unassigned))
		ZZ <- matrix(FALSE, n + n.groups, n - m)
		idx <- cumsum(sapply(groups, length)) + 1:n.groups
		ZZ[-idx, ] <- zz

		image(x = M,
					y = 1:(n + n.groups),
					z	= t(ZZ),
					zlim = c(0,n.groups + 1),
					col = plot.colors[1:(n.groups + 2)],
					xlab = x.label,
					ylab = y.label,
					main = main.label,
					axes = FALSE)

		ticks.at <- c(0.5,
		              cumsum(sapply(groups, length)) + 1:n.groups,
									n + n.groups + 0.5)
		labels.at <- (ticks.at[-1] + ticks.at[-length(ticks.at)]) / 2.0
		labels <- paste(group.names, ":  ", sapply(groups, length), " units", sep = "")
		labels <- c(labels, "unassigned")
		axis(2, at = ticks.at, labels = FALSE)
		axis(2, at = labels.at, labels = labels, tick = FALSE)

		abline(h = idx)
		axis(1)
		box()
	}
	
	else {

		image(x = M,
					y = 1:n,
					z	= t(zz),
					zlim = c(0,1),
					col = colors()[c(1,260)],
					xlab = x.label,
					ylab = y.label,
					main = main.label)
	
	  box()	
	}

  invisible(ans)
}

