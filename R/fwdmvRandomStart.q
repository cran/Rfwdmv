fwdmvRandomStart <- function(X, n.starts = 50, scaled = TRUE, plot.it = TRUE)
{
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  storage.mode(X) <- "double"

  starts <- matrix(0, n.starts, p+1)
  for(i in 1:n.starts)
    starts[i, ] <- sample(1:n, p+1)
  storage.mode(starts) <- "integer"

  mpomat <- matrix(as.double(0.0), n-p-1, n.starts)
  storage.mode(mpomat) <- "double"

  out <- .C("FSfwdmvRandomStart",
            X = X,
            n = as.integer(n),
            p = as.integer(p),
            starts = starts,
            n.starts = as.integer(n.starts),
            pbsb = as.integer(p+1),
            scaled = as.integer(scaled),
            mpomat = mpomat,
            PACKAGE = "Rfwdmv")

  if(n-p > 20) {
    domain <- (p+6):(n-1)
    out$mpomat <- out$mpomat[-(1:5), ]
  }
  else
    domain <- (p+1):(n-1)

  value <- list(starts = starts,
                mpomat = out$mpomat,
                domain = domain,
                n = n,
                p = p,
                m = min(domain),
                scaled = scaled)
 
  if(plot.it)
    plot.fwdmvRandomStart(value)

  oldClass(value) <- "fwdmvRandomStart"
  invisible(value)
}


