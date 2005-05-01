fwdmv <- function(X, groups = NULL, alpha = 0.6, beta = 0.75, bsb = ellipse.subset,
                  balanced = TRUE, scaled = TRUE, constrained = TRUE)
{
  data.name <- deparse(substitute(X))
  the.call <- match.call()

  if(class(X) == "fwdmv") {
    groups <- X$groups
    n.groups <- length(groups)
    group.names <- X$group.names
    x.names <- X$data.names
    X <- X$data
    dimnames(X) <- NULL
  }
  
  else {
    if(is.null(groups)) {
      ans <- fwdmv.init(X, bsb = bsb, scaled = scaled)
      ans$call <- the.call
      ans$data.name <- data.name
      return(ans)
    }

    X <- as.matrix(X)
    x.names <- dimnames(X)
    dimnames(X) <- NULL
    n.groups <- length(groups)

    if(!is.null(names(groups)))
      group.names <- names(groups)
    else
      group.names <- paste("Group", 1:n.groups)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]
  
  n.groups <- length(groups)
  l.groups <- sapply(groups, length)
  constrain <- constrained

  n.bsb <- ceiling(alpha * l.groups)
  n.bsb <- sapply(n.bsb, function(u, p) max(c(u, p)), p = p)
  M <- sum(n.bsb)

  bsb.threshold <- ceiling(beta * l.groups)
  unassigned <- setdiff(1:n, unlist(groups))
  n.unassigned <- length(unassigned)

  Cov <- list()
  Determinant <- list()
  Unit <- list()
  Center <- list()

  Max <- rep(as.double(NA), n-M+1)
  Mth <- rep(as.double(NA), n-M+1)
  Min <- rep(as.double(NA), n-M)
  Mpo <- rep(as.double(NA), n-M)

  for(i in 1:n.groups) {
    Center[[i]] <- matrix(0.0, n-M+1, p)
    Cov[[i]] <- matrix(0.0, n-M+1, (p*(p+1))/2)
    Determinant[[i]] <- double(n-M+1)
  }

  Distances <- matrix(0.0, n, n-M+1)
  next.dist <- double(n.groups)
  names(next.dist) <- 1:n.groups

  distances <- matrix(0.0, n, n.groups)
  dimnames(distances) <- list(1:n, 1:n.groups)

  Nearest <- matrix(as.integer(0), n.unassigned, n-M+1)
  dimnames(Nearest) <- list(unassigned, M:n)
  
  Misclassified <- list()

  # Find initial subsets in each group. #

  if(is.function(bsb)) {

    tmp.fun <- bsb
    bsb <- list()

    for(i in 1:n.groups)
      bsb[[i]] <- groups[[i]][tmp.fun(X[groups[[i]], , drop = FALSE], n.bsb[i])]

  }

  # Main forward search loop. #

  for(m in M:n) {
  
    Unit[[m-M+1]] <- bsb

    # Estimate the center and compute the R-factor of the matrices #
    # corresponding to the units in each group.                    #

    for(i in 1:n.groups) {

      center <- apply(X[bsb[[i]], , drop = FALSE], 2, mean)
      Xg <- sweep(X, 2, center)
      R <- qr.R(qr(Xg[bsb[[i]], ]))

      # Compute the covariance matrix and determinant for the group. #

      determinant <- prod(diag(R))^2 / (n.bsb[[i]] - 1)^p
      cov <- (t(R) %*% R) / (n.bsb[[i]] - 1)

      # Compute Mahalanobis distances for the units in each group. #

      u <- forwardsolve(t(R), t(Xg))
      u <- sqrt((n.bsb[[i]]-1) * colSums(u * u))

      if(scaled)
        distances[, i] <- u * determinant^(1/(2*p))
      else
        distances[, i] <- u

      if(!balanced)
        next.dist[i] <- ifelse(n.bsb[i]+1 <= l.groups[i], sort(distances[groups[[i]], i])[n.bsb[i]+1], Inf)

      # Store the center, covariance matrix (in compact form) and #
      # determinant for each group at this step.                  #

      Center[[i]][m-M+1, ] <- center
      Cov[[i]][m-M+1, ] <- cov[row(cov) <= col(cov)]
      Determinant[[i]][m-M+1] <- determinant

    }

    nearest <- apply(distances, 1, function(u) which(u == min(u)))
    Nearest[, m-M+1] <- nearest[unassigned]

    for(i in 1:n.groups) {
      if(any(misclassified <- nearest[groups[[i]]] != i)) {
        misclassified <- groups[[i]][misclassified]

        # store misclassified units: step, unit, nearest center #

        for(mis in misclassified) {
          mc <- m
					names(mc) <- nearest[mis]
          Misclassified[[as.character(mis)]] <- c(Misclassified[[as.character(mis)]], mc)
        }
      }
    }

    for(i in 1:n.groups)
      distances[groups[[i]], 1] <- distances[groups[[i]], i]

    distances[unassigned, 1] <- apply(distances[unassigned, , drop = FALSE], 1, min)

    sorted.distances <- sort(distances[, 1])

    Max[m-M+1] <- max(distances[unlist(bsb), 1])
    Mth[m-M+1] <- sorted.distances[m]



    if(m < n) {
      Min[m-M+1] <- min(distances[setdiff(1:n, unlist(bsb)), 1])
      Mpo[m-M+1] <- sorted.distances[m+1]
    }

    # Store the distances for all units. #

    Distances[, m-M+1] <- distances[, 1]

    # Update the subsets. #

    if(m < n) {

    # Start with constrained updating. #

      if(constrain) {
        if(balanced) {
          sizes <- n.bsb / l.groups
          i <- which(sizes == min(sizes))[1]
          n.bsb[i] <- n.bsb[i] + 1
          for(i in 1:n.groups)
            bsb[[i]] <- as.integer(names(sort(distances[groups[[i]], 1])[1:n.bsb[i]]))
        }

        if(!balanced) {
          not.full <- next.dist[n.bsb < bsb.threshold]
          i <- as.integer(names(which(not.full == min(not.full))))
          n.bsb[i] <- n.bsb[i] + 1
          for(i in 1:n.groups)
            bsb[[i]] <- sort(as.integer(names(sort(distances[groups[[i]], 1])[1:n.bsb[i]])))
        }

        if(all(n.bsb >= bsb.threshold))
          constrain <- FALSE
      }

    # Unconstrained updating. #

      else {

        # compute the subset #

        new.bsb <- as.integer(names(sort(distances[, 1])[1:(sum(n.bsb)+1)]))

        # compute within group subsets #

        for(i in 1:n.groups)
          bsb[[i]] <- sort(intersect(groups[[i]], new.bsb))

        # which unassigned units are in the subset #

        unassigned.in.subset <- intersect(unassigned, new.bsb)
        for(i in unassigned.in.subset)
          bsb[[nearest[i]]] <- sort(c(i, bsb[[nearest[i]]]))

        n.bsb <- sapply(bsb,length)
      }
    }
  }

  ans <- list(call = the.call,
              Distances = Distances,
              Center = Center,
              Cov = Cov,
              Determinant = Determinant,
              Unit = Unit,
              Nearest = Nearest,
              Misclassified = Misclassified,
              groups = groups,
              n = n,
              p = p,
              m = M,
              data = X,
              data.name = data.name,
              data.names = x.names,
              group.names = group.names,
              unassigned = unassigned,
              scaled = scaled,
              constrained = ifelse(constrained, as.integer(beta * sum(l.groups)), FALSE),
              Max = Max,
              Mth = Mth,
              Min = Min,
              Mpo = Mpo,
              initial = FALSE)

  oldClass(ans) <- "fwdmv"
  ans
}


