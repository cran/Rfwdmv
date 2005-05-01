bigunit.fwdmv <- function(x)
{
  n <- x$n
  m <- x$m

  bigunit <- matrix(FALSE, n, n-m+1)

  for(i in 1:(n-m+1))
    bigunit[unlist(x$Unit[[i]]), i] <- TRUE

  bigunit
}


