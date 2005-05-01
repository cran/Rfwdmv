fwdmvDistancePlot <- function(x, group = NULL, id = FALSE, psfrag.labels = FALSE)
{
  m <- x$m
  n <- x$n
  M <- m:n
  groups <- x$groups
	distances <- x$Distances

  ans <- list()
  color <- 16
  id.specified <- FALSE

  if(is.character(id) && casefold(id) == "all") {
		if(!is.null(group))
		  id <- groups[[group]]
		else
			id <- 1:n
  }

  if(is.numeric(id)) {
    selected <- unique(id)
    id.specified <- TRUE
    id <- TRUE
  }

  if(id && !id.specified && !interactive())
	  return(ans)

  if(id && !id.specified) {
    cat("\nLeft click to highlight trajectories. Right click\n")
    cat("to return to the R command prompt.\n\n")
  }

  if(!is.null(group)) {
    row.names <- groups[[group]]
    distances <- t(distances[row.names, , drop = FALSE])
  }

  else {
    row.names <- 1:n
    distances <- t(distances)
  }

  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }

  else {
    y.label <- ifelse(x$scaled, "Scaled Mahalanobis Distances", "Mahalanobis Distances")
    x.label <- "Subset Size"
    main.label <- ""
  }

  matplot(M, distances, type = "l", xlab = x.label,
    ylab = y.label, col = color, lty = 1, main = main.label)

  if(id) {

    if(id.specified) {
      for(sel in selected) {
        sel.idx <- which(row.names == sel)
        lines(M, distances[, sel.idx], col = 1, lwd = 2)
        text(n, distances[n-m+1, sel.idx], paste(" ", sel, sep = ""), adj = 0)
        text(m, distances[1, sel.idx], paste(sel, " ", sep = ""), adj = 1)
      }
    }

    else {
      counter <- 1
      selected <- character(0)
      while(counter < 20) {
        xy <- locator(1)
        if(is.null(xy))
          break

        xy$x <- round(xy$x) - m + 1
        step.dist <- abs(distances[xy$x, ] - xy$y)
        sel <- which(step.dist == min(step.dist))
        selected <- c(selected, row.names[sel])
        lines(M, distances[, sel], col = 1, lwd = 2)
        text(n, distances[n-m+1, sel], paste(" ", row.names[sel], sep = ""), adj = 0)
        text(m, distances[1, sel], paste(row.names[sel], " ", sep = ""), adj = 1)

        counter <- counter + 1
      }
    }

    ans$selected <- as.integer(unique(selected))
  }

  invisible(ans)
}

