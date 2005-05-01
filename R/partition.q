partition <- function(x, group = "next")
{
  m <- x$m
  n <- x$n
  p <- x$p
  M <- m:n
  unassigned <- x$unassigned
  groups <- x$groups
  n.groups <- length(groups)
  group.names <- x$group.names
	data.names <- x$data.names
  distances <- x$Distances
	dimnames(distances) <- list(data.names[[1]], M)

  if(!interactive())
	  return(x)

  max.dist <- max(distances)
  u.distances <- distances[unassigned, , drop = FALSE]

  y.label <- ifelse(x$scaled, "Scaled Mahalanobis Distances", "Mahalanobis Distances")
  x.label <- "Subset Size"
  main.label <- ""

  plot(M,
       rep(0.0, n-m+1),
       type = "n",
       xlab = x.label,
       ylab = y.label,
       main = main.label,
       ylim = c(0.0, max.dist))

  if(length(u.distances))
    matlines(M, t(u.distances), col = 8, lty = 1)

  if(n.groups) {
    for(i in 1:n.groups) {
      group.median <- apply(distances[groups[[i]], , drop = FALSE], 2, median)
      lines(M, group.median, lty = i, lwd = 2)
    }

    labels <- paste(group.names, "(", sapply(groups, length), " units)", sep = "")

    legend(n,
           max.dist,
           legend = labels,
           lty = 1:n.groups,
           lwd = 2,
           xjust = 1)
  }

  xy1 <- locator(1)
  xy1$x <- round(xy1$x)
  points(xy1$x, xy1$y, pch = 16, col = 2)

  xy2 <- locator(1)
  xy2$x <- round(xy2$x)
  points(xy2$x, xy2$y, pch = 16, col = 2)
  lines(c(xy1$x, xy2$x), c(xy1$y, xy2$y), lwd = 2, col = 2)
  
  pts <- c(xy1, xy2)

  xy1$x <- xy1$x - m + 1
  xy2$x <- xy2$x - m + 1


	if(xy2$x == xy1$x)
	  slope <- 1.0
	else
		slope <- sign((xy2$y - xy1$y) / (xy2$x - xy1$x))

	if(slope >= 0.0) {

	  if(xy1$y > xy2$y) {
		  temp <- xy2
			xy2 <- xy1
			xy1 <- temp
		}

		units <- as.integer(names(which(u.distances[, xy1$x] > xy1$y & u.distances[, xy2$x] < xy2$y)))
	}
	
	else {
	
	  if(xy1$y < xy2$y) {
		  temp <- xy2
			xy2 <- xy1
			xy1 <- temp
		}
	
		units <- as.integer(names(which(u.distances[, xy1$x] < xy1$y & u.distances[, xy2$x] > xy2$y)))
	}

	if(group == "next")
		group <- length(x$groups) + 1

  if(group == length(x$groups) + 1) {
    x$groups <- c(x$groups, list(units))
    x$group.names[group] <- paste("Group", group)
  }

  else {
    x$groups[[group]] <- c(x$groups[[group]], units)
  }

  x$unassigned <- setdiff(1:n, unlist(x$groups))
  x$pts <- pts
  x$initial <- FALSE

	fwdmvPartitionPlot(x, pts)

  x
}


