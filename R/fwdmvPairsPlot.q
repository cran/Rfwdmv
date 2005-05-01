fwdmvPairsPlot <- function(x)
{
  n <- x$n
  p <- x$p

  groups <- x$groups
  unassigned <- x$unassigned

  if(!length(groups)) {
    groups <- list(1:n)
    unassigned <- integer(0)
  }
  
  n.groups <- length(groups)
  n.in.groups <- sapply(groups, length)
  groups <- unlist(groups)

  ans <-list()

  if(n.groups <= 3)
    pch <- c(3,16,2)[1:n.groups]
  else
    pch <- 1:length(x$groups)

  pch <- rep(pch, n.in.groups)

  old.par <- par(mfrow = c(p, p), pty = "s", mar = rep(0.25,4), oma = rep(4, 4))
  on.exit(par(old.par))

  ranges <- apply(x$data, 2, range)

  for(i in 1:p) {
    for(j in 1:p) {
      if(i == j) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = ranges[,i], ylim = ranges[,i])
        text(mean(ranges[,j]), mean(ranges[,i]), x$data.names[[2]][i], cex = 2)
      }
      else {
        plot(x$data[groups, c(j, i)], pch = pch, axes = FALSE, xlim = ranges[,j], ylim = ranges[,i])
        if(length(unassigned))
          text(x$data[unassigned, c(j, i)], labels = unassigned, cex = 1.5)
      }

      box()
      if(i == p && j%%2) 
        axis(1)
      if(j == 1 && !(i%%2)) 
        axis(2)
      if(i == 1 && !(j%%2)) 
        axis(3)
      if(j == p && i%%2) 
        axis(4)
    }
  }
  
  invisible(ans)
}









