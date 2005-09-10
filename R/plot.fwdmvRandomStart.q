plot.fwdmvRandomStart <- function(x, envelope = TRUE, ...)
{
  n <- x$n
  p <- x$p
  m <- x$m
  domain <- x$domain
  mpomat <- x$mpomat
  scaled <- x$scaled

  if(envelope) {
    env <- FScomputeEnvelope(n, p, m)
    correction <- n*p*(domain-2)/((n-1)*(domain-p-1))
    env <- sqrt(sweep(env, 1, correction, "*"))

    if(!scaled) {
      c1 <- (domain-1)/n
      c2 <- qchisq(c1, p)
      correction <- sqrt(pchisq(c2, p+2)/c1)
      env <- sweep(env, 1, correction, "/")
    }

    y.range <- range(c(mpomat, env))
  }
  else
    y.range <- range(mpomat)

  matplot(domain,
          mpomat,
          ylim = y.range,
          type = "l",
          col = "gray",
          xlab = "Subset Size",
          ylab = "")

  if(envelope)
    matlines(domain, env, lty = 2, col = 1)

  invisible(x)
}
