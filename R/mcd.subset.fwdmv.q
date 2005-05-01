mcd.subset <- function(X, size)
{
  if(!is.element("package:lqs", search()))
    library(lqs)

  Xmcd <- cov.mcd(X)
  dist <- mahalanobis(X, Xmcd$center, Xmcd$cov)
  names(dist) <- 1:dim(X)[1]
  dist <- sort(dist)
  as.integer(names(dist[1:size]))
}


