panel.me <- function(x, y, scale = 1.0)
{
  medcov <- function(u, v) {
    n <- length(u)
    u <- u - median(u)
    v <- v - median(v)
    sum(u * v) / (n - 1)
  }
  
  sd <- sqrt(c(medcov(x, x), medcov(y, y)))
  phase <- acos(medcov(x, y) / (prod(sd)))
  center <- c(median(x), median(y))
  theta <- seq(-pi, pi, len = 101)
  xs <- center[1] + sd[1] * scale * cos(theta)
  ys <- center[2] + sd[2] * scale * cos(theta + phase)
  list(x = xs, y = ys)
}


