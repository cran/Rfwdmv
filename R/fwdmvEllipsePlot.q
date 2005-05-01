fwdmvEllipsePlot <- function(x, subset.size, plot.diagonal = TRUE)
{
  n <- x$n
  p <- x$p
  m <- x$m
  data <- x$data
  Distances <- x$Distances
  is.initial <- x$initial

  M <- m:n

	if(missing(subset.size))
		subset.size <- n

  ans <- list()

  step <- which(subset.size == M)

  groups <- x$groups
  unassigned <- x$unassigned

  if(is.initial) {
    groups <- x$Unit[[step]]
    unassigned <- integer(0)
  }

  group.c <- list()
  n.groups <- length(groups)

  if(n.groups <= 3)
    pch <- c(3,16,2)[1:n.groups]
  else
    pch <- 1:n.groups

  compute.ellipse <- function(center, cov, scale)
  {
    sd <- sqrt(diag(cov))
    phase <- acos(cov[1, 2]/(prod(sd)))
    theta <- seq(-pi, pi, len = 101)
    xs <- center[1] + sd[1] * scale * cos(theta)
    ys <- center[2] + sd[2] * scale * cos(theta + phase)
    list(x = xs, y = ys)
  }

  centers <- list()
  covs <- list()
  distances <- list()
  step.subset <- x$Unit[[step]]

  for(k in 1:n.groups) {
    centers[[k]] <- x$Center[[k]][step, ]

    covs[[k]] <- matrix(0.0, p, p)
    covs[[k]][row(covs[[k]]) <= col(covs[[k]])] <- x$Cov[[k]][step, ]
    covs[[k]][row(covs[[k]]) > col(covs[[k]])] <- t(covs[[k]])[row(covs[[k]]) > col(covs[[k]])]

    distances[[k]] <- Distances[step.subset[[k]], step]

    group.c[[k]] <- setdiff(groups[[k]], step.subset[[k]])
  }

  if(is.initial)
    group.c[[1]] <- setdiff(1:n, groups[[1]])

  pch <- rep(pch, sapply(group.c, length))
  group.c <- unlist(group.c)

  ### Compute Ellipses ###

  ellipses <- list()
  ec <- 1
  ranges <- apply(data, 2, range)

  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      for(k in 1:n.groups) {
        factor <- sqrt(max(distances[[k]]))
        ellipses[[ec]] <- compute.ellipse(centers[[k]][c(j,i)], covs[[k]][c(j,i), c(j,i)], factor)
        ranges[,i] <- range(c(ranges[,i], ellipses[[ec]]$y))
        ranges[,j] <- range(c(ranges[,j], ellipses[[ec]]$x))
        ec <- ec + 1

        factor <- sqrt(median(distances[[k]]))
        ellipses[[ec]] <- compute.ellipse(centers[[k]][c(j,i)], covs[[k]][c(j,i), c(j,i)], factor)
        ec <- ec + 1
      }
    }
  }

  ### Plot Ellipses ###
  
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

        bpdata <- list()
        for(k in 1:n.groups)
          bpdata[[k]] <- data[groups[[k]], i]
        boxplot(bpdata, axes = FALSE, ylim = ranges[,i])
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

        plot(data[group.c, j],
             data[group.c, i],
             xlim = ranges[,j],
             ylim = ranges[,i],
             pch = pch,
             xlab = "",
             ylab = "",
             axes = FALSE,
             type = "p")

        for(k in 1:n.groups) {
          lines(ellipses[[ec]], lty = k)
          ec <- ec + 1
          lines(ellipses[[ec]], lty = k)
          ec <- ec + 1
        }

        step.unassigned <- setdiff(unassigned, unlist(x$Unit[[step]]))
        if(length(step.unassigned))
          text(data[step.unassigned, c(j, i)], labels = step.unassigned, cex = 1.5)

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

  invisible(ans)
}



