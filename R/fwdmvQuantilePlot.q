fwdmvQuantilePlot <- function(x, subset.size, probs = "default", page = 1)
{
  if(probs[1] == "default")
    probs <- c(0.025, 0.05, 0.125, 0.25, 0.50, 0.75, 0.875, 0.95, 0.975)

  n <- x$n
  p <- x$p
  m <- x$m
  constrained <- x$constrained
  Distances <- x$Distances


  ans <- list()

  if(constrained < m)
    constrained <- FALSE

  M <- m:n
  step <- which(M == subset.size)

  distances <- Distances[unlist(x$Unit[[step]]), ]

  quan <- apply(distances, 2, quantile, probs = probs)

  step <- step + (page-1)*5
  if(step+5 > n-m+1)
    step <- n-m-4

  old.par <- par(mfrow = c(2,3), pty = "s")
  on.exit(par(old.par))

  units <- list()
  for(i in 1:5)
    units[[i]] <- setdiff(unlist(x$Unit[[step+i]]), unlist(x$Unit[[step+i-1]]))
  
  y.limit <- c(0.0, max(c(Distances[unlist(units), ]), distances))

  matplot(M,
          t(distances),
					type = "l",
					xlab = paste("Subset Size", M[step]),
					ylab = "",
					ylim = y.limit,
					main = "Reference Trajectories",
					col = 16)

  if(constrained)
    abline(v = constrained, col = 16, lty = 2)

  k <- 1
  for(i in (step+1):(step+5)) {

    unit <- units[[k]]
    k <- k + 1
    
    matplot(M,
		        t(quan),
		        type = "l",
						col = 16,
						xlab = paste("Subset Size", M[i]),
						ylab = "",
						main = paste(paste(unit, collapse = ", "), "Included"),
						ylim = y.limit)

    for(j in 1:length(unit))
      lines(M, Distances[unit[j], ])

    if(constrained)
      abline(v = constrained, col = 16, lty = 2)
  }

  invisible(ans)
}

