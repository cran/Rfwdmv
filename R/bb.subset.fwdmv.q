bb.subset <- function(X, size) {


bsplinex <- function(x, y)
{
# B-Spline helper function

	xr <- c(x, x[1:3])
	yr <- c(y, y[1:3])
	X <- rep(0, 10000)
	Y <- rep(0, 10000)
	zz <- 0
	storage.mode(xr) <- "double"
	storage.mode(yr) <- "double"
	storage.mode(X) <- "double"
	storage.mode(Y) <- "double"
	storage.mode(zz) <- "integer"
	xall <- .C("b_spline",
		x = xr,
		X = X,
		y = yr,
		Y = Y,
		z = zz,
		n = length(xr),
		PACKAGE = "Rfwdmv")
	splx <- xall$X[2:(xall$z - 30)]
	sply <- xall$Y[2:(xall$z - 30)]
	list(x=splx, y=sply)
}


outliers <- function(x, y, wx, wy, mx, my, eps)
{
# Outliers helper function
	
	X <- rep(0, length(x))
	Out <- rep(0, length(x))
	mm<-0
	nn<-length(x)
	ll = length(wx)
	storage.mode(X) <- "double"
	storage.mode(Out) <- "double"
	storage.mode(eps) <- "double"	
	storage.mode(x) <- "double"
	storage.mode(y) <- "double"
	storage.mode(wx) <- "double"
	storage.mode(wy) <- "double"
	storage.mode(mx) <- "double"
	storage.mode(my) <- "double"
	storage.mode(nn) <- "integer"
        storage.mode(ll) <- "integer"
	storage.mode(mm) <- "integer"
	soluz <- .C("outcor",
		X = X,
		Out = Out,
		eps = eps,
		x = x,
		y = y,
		wx = wx,
		wy = wy,
		mx = mx,
		my = my,
		nn = length(x),
		ll = length(wx),
		mm = mm,
		PACKAGE = "Rfwdmv")
	outlx <- soluz$X	
	ox <- which(soluz$X != 0)
	list(x=outlx, px=ox)
}


fence <- function(splx,sply,medx,medy,rad){
# Fence helper function

	XG <- (splx - medx)
	YG <- (sply - medy)

	XGK <- XG * rad
	YGK <- YG * rad
	
	xb <- (XGK + medx)
	yb <- (YGK + medy)
	
	list(x = xb, y = yb)
}




hull <- function(x, r = 1.0) {
# Hull helper function

	X <- as.matrix(x)
	XX <- X
	PX <- NULL
	
       
	box <- round(dim(X)[1] / 2)
	while (dim(X)[1] >= box) {
	     		hpts <- chull(X)
	     		X <- X[-hpts, ]
	     }
     
	hpts <- chull(X)
	hpts <- c(hpts, hpts[1])
	spl <- bsplinex(X[hpts, 1], X[hpts, 2])

	medx <- median(X[hpts, 1])
	medy <- median(X[hpts, 2])

	fen <- fence(spl$x, spl$y, medx, medy, r)
	outl <- outliers(XX[,1], XX[,2], fen$x, fen$y, medx, medy, 0.001)
	
	px <- outl$px

}

# main()

r <- 4.0


n <- nrow(X)
p <- ncol(X)

#Fr <- matrix(0, n, 2)

Fr <- matrix(0, n, 3)

Fr[ ,1] <- seq(1, n)
Fr[ ,3] <- 1

tmp <- Fr

if (p == 2) {
	while(sum(Fr[, 3]) > size) {
	Fr <- tmp
		outl <- hull(X, r)
		Fr[outl, 2] <- Fr[outl, 2] + 1
		Fr[outl, 3] <- 0
	r <- r - 0.1
}	
	
	sFr <- order(Fr[,2])
	bsb <- Fr[sFr[1:size], 1]
	
	bsb
	

} 

else if(p > 2) { 

   while(sum(Fr[, 3]) > size) {
   
	bsb <- seq(1, n)
	Fr <- tmp
	
	for (i in 1:(dim(X)[2] - 1)) {
		for (j in (i + 1):dim(X)[2]) {
			outl <- hull(cbind(X[,i], X[,j]), r)
			Fr[outl, 2] <- Fr[outl, 2] + 1
			Fr[outl, 3] <- 0
			#bsb <- setdiff(bsb, outl)
		
			}
	
	     }
	r <- r - 0.1
}

sFr <- order(Fr[,2])
bsb <- Fr[sFr[1:size], 1]
bsb
#Fr
  }

}

