fwdmv.init <- function(X, bsb = ellipse.subset, n.bsb, scaled = TRUE, monitor = "all")
{
  the.call <- match.call()
  data.name <- deparse(substitute(X))

  monitor.arg <- monitor
  monitor <- logical(9)
  names(monitor) <- c("distance", "center", "cov", "determinant",
                      "unit", "max", "mth", "min", "mpo")

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
  n <- dim(X)[1]
  p <- dim(X)[2]

  if(missing(n.bsb))
    n.bsb <- as.integer(p+1)

  if(is.function(bsb))
    bsb <- bsb(X, n.bsb)

  nbsb <- length(bsb)
  steps <- n-nbsb+1
  longbsb <- 1:n
  longbsb[1:nbsb] <- bsb

  storage.mode(X) <- "double"

  if(monitor["distance"]) {
    Distances <- matrix(as.double(0.0), n, steps)
    storage.mode(Distances) <- "double"
  }
  else
    Distances <- double(1)
    
  if(monitor["center"]) {
    Center <- matrix(as.double(0.0), steps, p)
    storage.mode(Center) <- "double"
  }
  else
    Center <- double(1)

  if(monitor["cov"]) {
    Cov <- matrix(as.double(0.0), steps, (p*(p+1))/2)
    storage.mode(Cov) <- "double"
  }
  else
    Cov <- double(1)

  if(monitor["determinant"])
    Determinant <- double(steps)
  else
    Determinant <- double(1)

  if(monitor["unit"]) {
    Unit <- matrix(as.integer(0), n, steps)
    storage.mode(Unit) <- "integer"
  }
  else
    Unit <- integer(1)

  if(monitor["max"])
    Max <- double(steps)
  else
    Max <- double(1)

  if(monitor["mth"])
    Mth <- double(steps)
  else
    Mth <- double(1)

  if(monitor["min"])
    Min <- double(steps - 1)
  else
    Min <- double(1)

  if(monitor["mpo"])
    Mpo <- double(steps - 1)
  else
    Mpo <- double(1)

  # working space for the native code #
  
  QR <- matrix(as.double(0.0), n, p)
  storage.mode(QR) <- "double"
  A <- matrix(as.double(0.0), n, p)
  storage.mode(A) <- "double"
  center <- double(p)
  distances <- double(n)
  wrk1n <- double(n)
  wrk2n <- double(n)
  newbsb <- integer(n)

  out <- .C("FSfwdmv",
            X = X,
            n = as.integer(n),
            p = as.integer(p),
            bsb = as.integer(longbsb - 1),
            n.bsb = as.integer(nbsb),
            monitor = as.integer(monitor),
            scaled = scaled,
            Distances = Distances,
            Center = Center,
            Cov = Cov,
            Determinant = Determinant,
            Unit = Unit,
            Max = Max,
            Mth = Mth,
            Min = Min,
            Mpo = Mpo,
            QR = QR,
            A = A,
            center = center,
            distances = distances,
            wrk1n = wrk1n,
            wrk2n = wrk2n,
            newbsb = newbsb,
            PACKAGE = "Rfwdmv")

  value <- list(call = the.call,
                Distances = if(monitor["distance"]) out$Distances else NA,
                Center = if(monitor["center"]) list(out$Center) else NA,
                Cov = if(monitor["cov"]) list(out$Cov) else NA,
                Determinant = if(monitor["determinant"]) list(out$Determinant) else NA,
                Unit = if(monitor["unit"]) out$Unit else NA,
                groups = list(),
                n = n,
                p = p,
                m = nbsb,
                data = X,
                data.name = data.name,
                data.names = x.names,
                group.names = "Unassigned",
                unassigned = 1:n,
                constrained = FALSE,
                scaled = scaled,
                Max = if(monitor["max"]) out$Max else NA,
                Mth = if(monitor["mth"]) out$Mth else NA,
                Min = if(monitor["min"]) out$Min else NA,
                Mpo = if(monitor["mpo"]) out$Mpo else NA,
                initial = TRUE)

  oldClass(value) <- "fwdmv"
  value
}



