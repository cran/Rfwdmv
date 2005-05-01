fwdmv.init <- function(X, bsb = ellipse.subset, scaled = TRUE)
{
  the.call <- match.call()
  data.name <- deparse(substitute(X))
  X <- as.matrix(X)
  x.names <- dimnames(X)
  dimnames(X) <- NULL

  n <- nrow(X)
  p <- ncol(X)

  n.bsb <- as.integer(nrow(X) / 3)

  if(is.function(bsb))
    bsb <- bsb(X, n.bsb)

  M <- length(bsb)

  Distances <- matrix(as.double(NA), nrow = n, ncol = n-M+1)
  Unit <- list()
  Center <- list(matrix(as.double(NA), nrow = n-M+1, ncol = p))
  Cov <- list(matrix(as.double(NA), nrow = n-M+1, ncol = (p*(p+1))/2))

  Max <- rep(as.double(NA), n-M+1)
  Mth <- rep(as.double(NA), n-M+1)
  Min <- rep(as.double(NA), n-M)
  Mpo <- rep(as.double(NA), n-M)
  Determinant <- list(rep(as.double(NA), n-M+1))

  newbsb <- integer(n)
  n1 <- 1:n

  for(m in M:n) {

    Unit[[m-M+1]] <- list(bsb)

    center <- apply(X[bsb, ], 2, mean)
    Center[[1]][m-M+1, ] <- center
    Xm <- sweep(X, 2, center)

    R <- qr.R(qr(Xm[bsb, ]))
    u <- forwardsolve(t(R), t(Xm))

    Xcov <- (t(R) %*% R)/(m-1)
    Cov[[1]][m-M+1, ] <- Xcov[row(Xcov) <= col(Xcov)]

    determinant <- prod(diag(R))^2 / (length(bsb) - 1)^p

    if(scaled)
      distances <- sqrt((m-1) * colSums(u * u)) * determinant^(1/(2*p))
    else
      distances <- sqrt((m-1) * colSums(u * u))

    Distances[, m-M+1] <- distances
    names(distances) <- n1
    Determinant[[1]][m-M+1] <- determinant

    sorted.distances <- sort(distances)

    Max[m-M+1] <- max(distances[bsb])
    Mth[m-M+1] <- sorted.distances[m]

    if(m < n) {
      Min[m-M+1] <- min(distances[setdiff(n1, bsb)])
      Mpo[m-M+1] <- sorted.distances[m+1]

      newbsb <- as.integer(names(sorted.distances))[1:(m+1)]

      bsb <- newbsb
    }

  }

  value <- list(call = the.call,
                Distances = Distances,
                Center = Center,
                Cov = Cov,
                Determinant = Determinant,
                Unit = Unit,
                groups = list(),
                n = n,
                p = p,
                m = M,
                data = X,
                data.name = data.name,
                data.names = x.names,
                group.names = "Unassigned",
                unassigned = n1,
                constrained = FALSE,
                scaled = scaled,
                Max = Max,
                Mth = Mth,
                Min = Min,
                Mpo = Mpo,
                initial = TRUE)

  oldClass(value) <- "fwdmv"
  value
}


