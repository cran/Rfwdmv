FScomputeEnvelope <- function(n, p, m, probs = c(0.01, 0.5, 0.99),
                              tolerance = 1e-3)
{
  nprobs <- length(probs)
  quantiles <- matrix(0.0, nrow = n-m, ncol = nprobs)
  storage.mode(quantiles) <- "double"

  ans <- .C("FScomputeEnvelope",
            n = as.integer(n),
            p = as.integer(p),
            m = as.integer(m),
            probs = as.double(probs),
            nprobs = as.integer(nprobs),
            quantiles = quantiles,
            tol = as.double(tolerance),
            PACKAGE = "Rfwdmv")

  dimnames(ans$quantiles) <- list(m:(n-1), paste(100*probs, "%", sep = ""))
  ans$quantiles
}


