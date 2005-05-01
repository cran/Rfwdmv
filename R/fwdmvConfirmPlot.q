fwdmvConfirmPlot <- function(x, n.steps, psfrag.labels = FALSE)
{
  n <- x$n
  m <- x$m
  is.initial <- x$initial
  Nearest <- x$Nearest
  Misclassified <- x$Misclassified
  groups <- x$groups
  unassigned <- x$unassigned

  n.groups <- length(groups)
  n.unassigned <- length(unassigned)
  misclassified <- sort(as.integer(names(Misclassified)))
  n.misclassified <- length(misclassified)

  ans <- list()

  if(is.initial)
    stop("This plot is not available for initial fits: run fwdmv first.")

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- ""
    x.label <- "Subset Size"
    main.label <- ""
  }

  respective.groups <- integer(n.misclassified)

  if(n.misclassified) {
		for(i in 1:n.misclassified) {
			for(j in 1:n.groups) {
				if(is.element(misclassified[i], groups[[j]])) {
					respective.groups[i] <- j
					break
				}
			}
		}

		bottom <- matrix(rep(respective.groups, each = n.steps), nrow = n.misclassified, byrow = TRUE)
		dimnames(bottom) <- list(misclassified, (n-n.steps+1):n)
		
		for(i in 1:n.misclassified) {
			mis <- Misclassified[[as.character(misclassified[i])]]
			mis <- mis[mis >= n-n.steps+1]
			bottom[i, mis+n.steps-n] <- as.integer(names(mis))
		}
	}

  top <- Nearest[, (n-n.steps+1):n - m + 1, drop = FALSE]

  if(n.misclassified)
    both <- rbind(bottom, top)
  else
	  both <- top

  if(n.groups <= 3)
    pch <- c(3,16,2)[1:n.groups]
  else
    pch <- 1:n.groups

  old.par <- par(mfrow = c(1,1))
  on.exit(par(old.par))

  both <- pch[t(both)]

  if(n.misclassified) {
		plot(expand.grid((n-n.steps+1):n, c(1:n.misclassified, (n.misclassified+2):(n.misclassified+n.unassigned+1))),
				 pch = both,
				 ylab = y.label,
				 xlab = x.label,
				 main = main.label,
				 ylim = c(1 - n.groups, n.unassigned + n.misclassified+2),
				 axes = FALSE)

		box()
		axis(1)
		axis(2,
				 at = c(1:n.misclassified, (n.misclassified+2):(n.misclassified+n.unassigned+1)),
				labels = c(misclassified, unassigned),
				las = 1,
				adj = 1)
		abline(h = 0.5)
		abline(h = n.misclassified + 1.5)
		text(mean(par()$usr[1:2]), n.misclassified+1, "Misclassified Units")
		text(mean(par()$usr[1:2]), n.unassigned + n.misclassified+2, "Unassigned Units")
		legend(n-n.steps+1, 0, legend = x$group.names, pch = pch, xjust = 0, yjust = 1, bty = "n",
					 y.intersp = 1.25, ncol = n.groups)
  }
	
	else {
		plot(expand.grid((n-n.steps+1):n, 1:n.unassigned),
				 pch = both,
				 ylab = y.label,
				 xlab = x.label,
				 main = main.label,
				 ylim = c(1 - n.groups, n.unassigned),
				 axes = FALSE)

		box()
		axis(1)
		axis(2,
				at = 1:n.unassigned,
				labels = unassigned,
				las = 1,
				adj = 1)
		abline(h = 0.5)
		text(mean(par()$usr[1:2]), 0, "Unassigned Units")
		legend(n - n.steps + 1,
		       0,
					 legend = x$group.names,
					 pch = pch,
					 xjust = 0,
					 yjust = 1,
					 bty = "n",
					 y.intersp = 1.25,
					 ncol = n.groups)
	}

  invisible(ans)
}






