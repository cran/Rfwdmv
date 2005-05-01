panel.be <- function(x, y, scale = 1.0)
{
  sd <- sqrt(c(var(x), var(y)))
  phase <- acos(var(x, y) / (prod(sd)))
  center <- c(mean(x), mean(y))
  theta <- seq(-pi, pi, len = 101)
  xs <- center[1] + sd[1] * scale * cos(theta)
  ys <- center[2] + sd[2] * scale * cos(theta + phase)
  list(x = xs, y = ys)
}


