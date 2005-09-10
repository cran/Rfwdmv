fwdmv.init.pureR <- function(X, bsb = ellipse.subset, n.bsb, scaled = TRUE, monitor = "all")
{
  the.call <- match.call()
  data.name <- deparse(substitute(X))

  monitor.arg <- monitor
  monitor <- logical(9)
  names(monitor) <- c("cov", "determinant", "unit", "center", "max", "mth",
                      "min", "mpo", "distance")

  if(casefold(monitor.arg) == "all")
    monitor[1:9] <- TRUE

  else
    for(stat in monitor.arg) {
      stat <- casefold(stat)
      monitor[stat] <- TRUE
    }

  X <- as.matrix(X)
  x.names <- dimnames(X)
  dimnames(X) <- NULL

  n <- nrow(X)
  p <- ncol(X)

  if(missing(n.bsb))
    n.bsb <- as.integer(p+1)

  if(is.function(bsb))
    bsb <- bsb(X, n.bsb)

  M <- length(bsb)

  if(monitor["distance"])
    Distances <- matrix(as.double(NA), nrow = n, ncol = n-M+1)
  if(monitor["unit"])
    Unit <- matrix(0, nrow = n, ncol = n-M+1)
  if(monitor["center"])
    Center <- list(matrix(as.double(NA), nrow = n-M+1, ncol = p))
  if(monitor["cov"])
    Cov <- list(matrix(as.double(NA), nrow = n-M+1, ncol = (p*(p+1))/2))

  if(monitor["max"])
    Max <- rep(as.double(NA), n-M+1)
  if(monitor["mth"])
    Mth <- rep(as.double(NA), n-M+1)
  if(monitor["min"])
    Min <- rep(as.double(NA), n-M)
  if(monitor["mpo"])
    Mpo <- rep(as.double(NA), n-M)
  if(monitor["determinant"])
    Determinant <- list(rep(as.double(NA), n-M+1))

  newbsb <- integer(n)
  n1 <- 1:n

  for(m in M:n) {

    if(monitor["unit"])
      Unit[bsb, m-M+1] <- 1

    center <- apply(X[bsb, ], 2, mean)
    if(monitor["center"])
      Center[[1]][m-M+1, ] <- center
    Xm <- sweep(X, 2, center)

    R <- qr.R(qr(Xm[bsb, ]))
    u <- forwardsolve(t(R), t(Xm))

    Xcov <- (t(R) %*% R)/(m-1)
    if(monitor["cov"])
      Cov[[1]][m-M+1, ] <- Xcov[row(Xcov) <= col(Xcov)]

    if(monitor["determinant"] || scaled)
      determinant <- prod(diag(R))^2 / (length(bsb) - 1)^p

    if(scaled)
      distances <- sqrt((m-1) * colSums(u * u)) * determinant^(1/(2*p))
    else
      distances <- sqrt((m-1) * colSums(u * u))

    if(monitor["distance"])
      Distances[, m-M+1] <- distances
    names(distances) <- n1
    if(monitor["determinant"])
      Determinant[[1]][m-M+1] <- determinant

    sorted.distances <- sort(distances)

    if(monitor["max"])
      Max[m-M+1] <- max(distances[bsb])
    if(monitor["mth"])
      Mth[m-M+1] <- sorted.distances[m]

    if(m < n) {
      if(monitor["min"])
        Min[m-M+1] <- min(distances[setdiff(n1, bsb)])
      if(monitor["mpo"])
        Mpo[m-M+1] <- sorted.distances[m+1]

      newbsb <- as.integer(names(sorted.distances))[1:(m+1)]
      bsb <- newbsb
    }
  }

  value <- list(call = the.call,
                Distances = if(monitor["distance"]) Distances else NA,
                Center = if(monitor["center"]) Center else NA,
                Cov = if(monitor["cov"]) Cov else NA,
                Determinant = if(monitor["determinant"]) Determinant else NA,
                Unit = if(monitor["unit"]) Unit else NA,
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
                Max = if(monitor["max"]) Max else NA,
                Mth = if(monitor["mth"]) Mth else NA,
                Min = if(monitor["min"]) Min else NA,
                Mpo = if(monitor["mpo"]) Mpo else NA,
                initial = TRUE)

  oldClass(value) <- "fwdmv"
  value
}


