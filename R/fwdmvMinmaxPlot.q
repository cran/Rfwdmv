fwdmvMinmaxPlot <- function(x, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  M <- m:n

  ans <-list()

  if(psfrag.labels) {
    x.label1 <- "xlab1"
    y.label1 <- "ylab1"
    main.label1 <- "main1"
    x.label2 <- "xlab2"
    y.label2 <- "ylab2"
    main.label2 <- "main2"
  }

  else {
    x.label1 <- "Subset Size"
    y.label1 <- "Maximum Included Distance and mth Distance"
    main.label1 <- ""
    x.label2 <- "Subset Size"
    y.label2 <- "Minimum Excluded Distance and (m+1)th Distance"
    main.label2 <- ""
  }

  old.par <- par(mfrow = c(1,2), pty = "s")
  on.exit(par(old.par))

  y.lim <- range(c(x$Max, x$Mth))

  plot(M,
       x$Max,
       type = "l",
       xlab = x.label1,
       ylim = y.lim,
       ylab = y.label1,
       main = main.label1,
			 lty = 2)

  lines(M, x$Mth)

  y.lim <- range(c(x$Min, x$Mpo))

  plot((m+1):n,
       x$Min,
       type = "l",
       xlab = x.label2,
       ylim = y.lim,
       ylab = y.label2,
       main = main.label2,
			 lty = 2)

  lines((m+1):n, x$Mpo)

  invisible(ans)
}


