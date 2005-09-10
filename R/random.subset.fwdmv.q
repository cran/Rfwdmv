random.subset <- function(X, n.bsb)
{
  n <- dim(X)[1]
  sample(1:n, n.bsb)
}

