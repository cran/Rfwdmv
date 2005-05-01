panel.bb <- function(x, y, scale = 1.0)
{
	bsplinex <- function(x, y) {

		x <- as.double(c(x, x[1:3]))
		y <- as.double(c(y, y[1:3]))
		X <- double(10000)
		Y <- double(10000)
		z <- as.integer(0)
		
		n <- as.integer(length(x))

		xall <- .C("b_spline",
			x = x,
			X = X,
			y = y,
			Y = Y,
			z = z,
			n = n,
			PACKAGE = "Rfwdmv")

		splx <- xall$X[2:(xall$z - 30)]
		sply <- xall$Y[2:(xall$z - 30)]
		list(x = splx, y = sply)
	}

  X <- cbind(x, y)
  X <- as.matrix(X)

  box <- round(dim(X)[1]/2)
     
  while(dim(X)[1] >= box) {
    hpts <- chull(X)
    X <- X[-hpts, ]
  }

  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])

  spl <- bsplinex(X[hpts, 1], X[hpts, 2])

  fence <- function(splx, sply, medx, medy, rad)
  {
    x <- (splx - medx) * rad + medx
    y <- (sply - medy) * rad + medy

    list(x = x, y = y)
  }

  medx <- median(X[hpts, 1])
  medy <- median(X[hpts, 2])

  if(scale != 1.0)
    spl <- fence(spl$x, spl$y, medx, medy, scale)

  spl
}


