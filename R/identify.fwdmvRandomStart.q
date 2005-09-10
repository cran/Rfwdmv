identify.fwdmvRandomStart <- function(x, ...)
{
  xy <- locator(1)
  if(is.null(xy))
    return(invisible(x))

  col.idx <- round(xy$x) - min(x$plot.domain) + 1
  nearest <- abs(x$mpomat[col.idx, ] - xy$y)
  trajectories <- which(nearest == min(nearest))
  for(trajectory in trajectories)
    lines(x$plot.domain, x$mpomat[, trajectory], col = 1, lwd = 2)

  value <- x$starts[trajectories, , drop = FALSE]
  dimnames(value) <- list(paste("start ", format(trajectories), ":", sep = ""), 1:dim(value)[2])
  value
}


