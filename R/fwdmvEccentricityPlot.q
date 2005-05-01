fwdmvEccentricityPlot <- function(x, which = c(1,2), psfrag.labels = FALSE)
{
  which <- sort(which)
  n <- x$n
  m <- x$m
  M <- m:n

  ans <- list()

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
    y.label1 <- "Fraction of Variance Explained"
    main.label1 <- ""

    x.label2 <- "Subset Size"
    y.label2 <- "Ecentricity"
    main.label2 <- ""
  }

  eigenvalues <- eigenvalues.fwdmv(x)
  n.groups <- length(eigenvalues)

  total.variance <- list()

  for(i in 1:n.groups) {
    total.variance[[i]] <- apply(eigenvalues[[i]], 1, sum)
    eigenvalues[[i]] <- eigenvalues[[i]][, which]
  }

  old.par <- par(mfrow = c(1,2), pty = "s")
  on.exit(par(old.par))

  plot(NA, NA,
       xlim = c(m, n),
       ylim = c(0,1),
       xlab = x.label1,
       ylab = y.label1,
       main = main.label1)

  for(i in 1:n.groups) {
    lines(M, eigenvalues[[i]][, 1] / total.variance[[i]], col = i)
    lines(M, eigenvalues[[i]][, 2] / total.variance[[i]], col = i, lty = 2)
  }

  plot(NA, NA,
       xlim = c(m, n),
       ylim = c(0,1),
       xlab = x.label2,
       ylab = y.label2,
       main = main.label2)

  for(i in 1:n.groups) {
    lines(M, 1 - sqrt(eigenvalues[[i]][,2]/eigenvalues[[i]][,1]), col = i)
  }

  invisible(ans)
}


