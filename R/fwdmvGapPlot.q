fwdmvGapPlot <- function(x, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
 
  ans <- list()

  if(psfrag.labels) {
    x.label1 <- "xlab1"
    y.label1 <- "ylab1"
    main.label1 <- "main1"
  }

  else {
    x.label1 <- "Subset Size"
    y.label1 <- "Difference"
    main.label1 <- "Gap Plot"
  }

  old.par <- par(pty = "s")
  on.exit(par(old.par))

  y1 <- x$Min - x$Max[-length(x$Max)]
  y2 <- x$Mpo - x$Mth[-length(x$Mth)]
  
  y.limits <- range(c(y1, y2))

  plot((m+1):n,
       y2,
       ylim = y.limits,
       type = "l",
       lwd = 4,
       xlab = x.label1,
       ylab = y.label1,
       main = main.label1,
			 col = 8)

  lines((m+1):n, y1, lty = 4)

  legend(m+1,
         par("usr")[4],
         legend = c("(m+1)th distance minus mth distance", "Min excluded distance minus max included distance"),
         col = c(8,1),
         lwd = c(4,2),
         lty = c(1,4))

  invisible(ans)
}


