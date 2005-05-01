eigenvectors.fwdmv <- function(x, which.vector = 1)
{
  p <- x$p
  n <- x$n
  m <- x$m

  Covariance <- x$Cov
  n.groups <- length(Covariance)
  eigenvectors <- list()

  tmp.fun1 <- function(u, p, which.vector) {
    mat <- matrix(0.0, p, p)
    mat[row(mat) <= col(mat)] <- u
    mat[col(mat) < row(mat)] <- t(mat)[col(mat) < row(mat)]
    eigen(mat)$vectors[, which.vector]
  }

  tmp.fun2 <- function(u, ms) {
    if(sign(u[abs(u) == max(abs(u))]) == ms)
      u <- -u
    u
  }

  for(k in 1:n.groups) {

    eigenvectors[[k]] <- t(apply(Covariance[[k]], 1, tmp.fun1, p = p, which.vector = which.vector))

    ms <- sign(eigenvectors[[k]][which(abs(eigenvectors[[k]]) == max(abs(eigenvectors[[k]])))])[1]
    eigenvectors[[k]][1, ] <- tmp.fun2(eigenvectors[[k]][1, ], ms)

    for(i in 2:(n-m+1)) {
      idx <- which(abs(eigenvectors[[k]][i,]) == max(abs(eigenvectors[[k]][i,])))
      if(sign(eigenvectors[[k]][i, idx]) != sign(eigenvectors[[k]][i-1, idx]))
        eigenvectors[[k]][i,] <- -1.0 * eigenvectors[[k]][i,]
    }
  }

  eigenvectors
}

