fwdmvPartitionPlot <- function(x, pts = NULL, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  p <- x$p
  M <- m:n
	data.names <- x$data.names
  distances <- x$Distances
	dimnames(distances) <- list(data.names[[1]], M)
  group.names <-x$group.names

  ans <- list()

  unassigned <- x$unassigned
  groups <- x$groups
  n.groups <- length(x$groups)
  
  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- ifelse(x$scaled, "Scaled Mahalanobis Distances", "Mahalanobis Distances")
    x.label <- "Subset Size"
    main.label <- ""
  }

  max.dist <- max(distances)
  u.distances <- distances[unassigned, , drop = FALSE]

  plot(M,
       rep(0.0, n-m+1),
       type = "n",
       xlab = x.label,
       ylab = y.label,
       main = main.label,
       ylim = c(0.0, max.dist))

  if(length(u.distances))
    matlines(M, t(u.distances), col = 8, lty = 1)
  
  if(!is.null(pts)) {
    xy1 <- pts[1:2]
    xy2 <- pts[3:4]
    points(xy1$x, xy1$y, pch = 16, col = 2)
    points(xy2$x, xy2$y, pch = 16, col = 2)
    lines(c(xy1$x, xy2$x), c(xy1$y, xy2$y), lwd = 2, col = 2)
  }

  if(n.groups) {
    for(i in 1:n.groups) {
      group.median <- apply(distances[groups[[i]], , drop = FALSE], 2, median)
      lines(M, group.median, lty = i, lwd = 2)
    }

    labels <- paste(group.names, " (", sapply(groups, length), " units)", sep = "")

    legend(n,
           max.dist,
           legend = labels,
           lty = 1:n.groups,
           lwd = 2,
           xjust = 1)
  }

  invisible(ans)
}


