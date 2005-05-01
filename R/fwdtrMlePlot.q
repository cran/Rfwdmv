fwdtrMlePlot <- function(x, psfrag.labels = FALSE)
{
  n <- x$n
  m <- x$m
  H0 <- x$H0
  Mle <- x$Mle

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }
  else{
    x.label <- "Subset Size"
    y.label <- ifelse(x$forced.onepar, "MLE of common Transformation Parameter", "MLE of Transformation Parameters")
    main.label <- " "
  }
  
  n.lambda <- length(H0[, 1])

  matplot(m:n,
          Mle[[1]],
          type = "l",
          col = 1,
          ylim = range(Mle[[1]]),
          lty = 1:n.lambda,
          lwd = 1:n.lambda,
          xlab = x.label,
          ylab = y.label,
          main = main.label
          )
  
  abline(h = -1.0, col = 8)
  abline(h = 0.0, col = 8)
  abline(h = 1.0, col = 8)

  if(!x$forced.onepar){
    text(n, Mle[[1]][n-m+1,], paste("  ", H0[,2], sep = ""), adj = 0, cex = 1.5)
    text(m, Mle[[1]][1,], paste(H0[,2], " ", sep = ""), adj = 1, cex = 1.5)
  }
  invisible()
}


