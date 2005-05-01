eigenvalues.fwdmv <- function(x)
{
  p <- x$p
  Covariance <- x$Cov
  n.groups <- length(Covariance)
  eigenvalues <- list()

  tmp.fun <- function(u, p) {
    mat <- matrix(0.0, p, p)
    mat[row(mat) <= col(mat)] <- u
    mat[col(mat) < row(mat)] <- t(mat)[col(mat) < row(mat)]
    eigen(mat)$values
  }

  for(i in 1:n.groups)
    eigenvalues[[i]] <- t(apply(Covariance[[i]], 1, tmp.fun, p = p))

  eigenvalues
}


