ellipse.subset <- function(X, size)
{
  n <- nrow(X)
  p <- ncol(X)

  medians <- apply(X, 2, median)
  Xm <- sweep(X, 2, medians)

  Distances <- matrix(as.double(NA), nrow = n, ncol = (p*(p-1)/2))
  counter <- 1

  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      x <- Xm[ , c(i,j)]
      R <- qr.R(qr(x))
      u <- forwardsolve(t(R), t(x))
      Distances[, counter] <- sqrt((n-1) * apply(u, 2, function(x) t(x) %*% x))
      counter <- counter + 1
    }
  }

  distances <- apply(Distances, 1, max)

  threshold <- quantile(distances, probs = (size / n))
  Subset <- which(distances <= threshold)

  Subset
}



