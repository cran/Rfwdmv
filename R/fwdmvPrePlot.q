fwdmvPrePlot <- function(X, panel = panel.be, plot.diagonal = TRUE)
{
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]

  ### Compute Contours ###

  contours <- list()
  ec <- 1
  ranges <- apply(X, 2, range)

  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      contours[[ec]] <- panel(X[, j], X[, i], 1.8)
      ranges[,i] <- range(c(ranges[,i], contours[[ec]]$y))
      ranges[,j] <- range(c(ranges[,j], contours[[ec]]$x))
      ec <- ec + 1

      contours[[ec]] <- panel(X[, j], X[, i], 1.0)
      ec <- ec + 1
    }
  }

  ### Plot Contours ###
  
  if(plot.diagonal) {
    rows <- 1:p
    cols <- 1:p
    mfrow <- c(p, p)
  }

  else {
    rows <- 1:(p-1)
    cols <- 2:p
    mfrow <- c(p-1, p-1)
  }

  old.par <- par(mfrow = mfrow, pty = "s", mar = rep(0.25,4), oma = rep(4, 4))
  on.exit(par(old.par))

  ec <- 1
  for(i in rows) {
    for(j in cols) {
    
      if(i == j && plot.diagonal) {
        boxplot(X[, i], axes = FALSE, ylim = ranges[,i], pch = 16)
        box()
      }

      else if(i >= j) {

        plot(mean(ranges[,j]),
             mean(ranges[,i]),
             xlim = ranges[,j],
             ylim = ranges[,i],
             xlab = "",
             ylab = "",
             axes = FALSE,
             type = "n")
      }

      else {

        plot(X[, j],
             X[, i],
             xlim = ranges[,j],
             ylim = ranges[,i],
             pch = 16,
             xlab = "",
             ylab = "",
             axes = FALSE,
             type = "p",
             col = 8)

        lines(contours[[ec]], lwd = 2)
        ec <- ec + 1
        lines(contours[[ec]], lwd = 2)
        ec <- ec + 1

        box()
      }

      if(j >= i) {
        if(plot.diagonal) {
          if(i == p && j%%2) 
            axis(1)
          if(j == 1 && !(i%%2)) 
            axis(2)
          if(i == 1 && !(j%%2)) 
            axis(3)
          if(j == p && i%%2) 
            axis(4)
        }

        else {
          if(i == p-1 && j%%2) 
            axis(1)
          if(j == 1 && !(i%%2)) 
            axis(2)
          if(i == 1 && !(j%%2)) 
            axis(3)
          if(j == p && i%%2) 
            axis(4)
        }
      }
    }
  }

  invisible()
}

