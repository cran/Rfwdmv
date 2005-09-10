mcd.subset <- function(X, size)
{
  Xmcd <- cov.mcd(X)
  dist <- mahalanobis(X, Xmcd$center, Xmcd$cov)
  names(dist) <- 1:dim(X)[1]
  dist <- sort(dist)
  as.integer(names(dist[1:size]))
}


