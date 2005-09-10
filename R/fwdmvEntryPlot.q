fwdmvEntryPlot <- function(x, entry.order = "first", subset.size = -1, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
	groups <- x$groups
	unassigned <- x$unassigned
	group.names <- x$group.names

  M <- m:n
  n.groups <- length(groups)
	n.unassigned <- length(unassigned)
  ans <- list()
	plot.colors <- colors()[c(1, 260, 387, 490, 536, 640, 498, 552, 652)]

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    x.label <- "Subset Size"
    y.label <- switch(entry.order,
      "first" = "First Entry Order",
      "final" = "Final Entry Order",
      "natural" = "Natural Entry Order",
      "integer" = paste("Step", subset.size, "Entry Order"))
    main.label <- ""
      
  }

  entry.order <- switch(entry.order,

    "first" = {
      j <- 1
      units <- numeric(0)
      while(length(units) < n) {
        units <- unique(c(units, which(x$Unit[ , j] > 0)))
        j <- j + 1
      }
      units
    },

    "final" = {
      j <- 1
      units <- numeric(0)
      while(length(units) < n) {
        units <- unique(c(which(x$Unit[ , j] > 0), units))
        j <- j + 1
      }
      units
    },
    
    "natural" = {
      1:n
    },

    "integer" = {
      units <- x$Distances[ , subset.size - m + 1]
      names(units) <- 1:n
      as.integer(names(sort(units)))
    }
  )

  if(length(groups)) {

		for(i in 1:n.groups)
			groups[[i]] <- intersect(entry.order, groups[[i]])

		unassigned <- intersect(entry.order, unassigned)
			
		entry.order <- c(unlist(groups), unassigned)	

		counts <- sapply(groups, length)
		cscounts <- cumsum(counts)
		bigunit <- matrix(0, n + n.groups, n - m + 1)
		bigunit[-(cscounts + 1:n.groups), ] <- x$Unit[entry.order, , drop = FALSE]
		group.id <- rep(0:n.groups + 1, times = (c(counts, n.unassigned) + 1))
		group.id <- group.id[-length(group.id)]
		bigunit <- bigunit * group.id

		image(x = M,
					y = 1:(n + n.groups),
					z	= bigunit,
					zlim = c(0, n.groups + 1),
					col = plot.colors[1:(n.groups + 2)],
					xlab = x.label,
					ylab = y.label,
					main = main.label,
					axes = FALSE)

		abline(h = cscounts + 1:n.groups)
		axis(1)

		ticks.at <- c(0.5,
									cumsum(sapply(groups, length)) + 1:n.groups,
									n + n.groups + 0.5)
		labels.at <- (ticks.at[-1] + ticks.at[-length(ticks.at)]) / 2.0
		labels <- paste(group.names, ":  ", sapply(groups, length), " units", sep = "")
		labels <- c(labels, "unassigned")
		axis(2, at = ticks.at, labels = FALSE)
		axis(2, at = labels.at, labels = labels, tick = FALSE)

		box()
  }

  else {

		bigunit <- x$Unit[entry.order, , drop = FALSE]
		
		image(x = M,
					y = 1:n,
					z	= t(bigunit),
					zlim = c(0, 1),
					col = plot.colors[1:2],
					xlab = x.label,
					ylab = y.label,
					main = main.label)

		box()
	}

  invisible(ans)
}


