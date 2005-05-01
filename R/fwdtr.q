fwdtr <- function(X, bsb = ellipse.subset, n.bsb = 50 , lambda = 1,
                        one.lambda = FALSE, col.to.transform = "all",
                         boundaries = c(-3, 3))
{
  X <- as.matrix(X)
  x.names <- dimnames(X)
  dimnames(X) <- NULL
  n <- nrow(X)
  p <- ncol(X)
  the.call <- match.call()
  if(is.character(col.to.transform)) {
    col.to.transform <- seq(1:ncol(X))
    boundaries <- matrix(rep(boundaries, ncol(X)), ncol = 2, byrow = TRUE)
  }
  # Check initial conditions and consistency of user's choices
  else {
    col.to.transform <- as.integer(col.to.transform)
    boundaries <- matrix(rep(boundaries, length(col.to.transform)), ncol = 2, byrow = TRUE)
  }

  if(any(X[,col.to.transform] <= 0.0))
    stop("The Box-Cox transformation is only defined for positive data.")

  if(length(lambda) > length(col.to.transform))
    stop("You have to select a wider number of columns to transform.")

  if((length(lambda) > 1) && (length(lambda) < length(col.to.transform)))
    stop("You must input a value of lambda for each column you have chosen.")

  if((length(lambda) == 1) && (length(col.to.transform) > 1))
    lambda <- rep(lambda, length(col.to.transform))

  new.lambda <- cbind(lambda, col.to.transform, boundaries)
  n.bsb <- max((p+1), as.integer((n.bsb * nrow(X))/100))
   if(is.function(bsb))
    bsb <- bsb(X, n.bsb)
  j <- 0
  X.lambda <- X  
  # Transform selected columns (ONLY) following BOX-COX under H0.   #
  # X.lambda is used for sorting distances and build Clean Data Set #
  for(ii in new.lambda[,2]) {
    j <- j+1
    G <- exp(mean(log(X[,ii])))
    if(new.lambda[j,1] == 0) {
      X.lambda[, ii] <- log(X[, ii]) * G
    }
    else {
      X.lambda[,ii] <- ((X[,ii]^new.lambda[j, 1]) - 1)/(new.lambda[j, 1] * G^(new.lambda[j, 1]-1))
    }
  }


  M <- length(bsb)
  Unit <- list()
  if(one.lambda)
    Mle <- list(matrix(as.double(NA), nrow=n-M+1, ncol = 1))
  else
    Mle <- list(matrix(as.double(NA), nrow=n-M+1, ncol = length(new.lambda[, 1])))

  newbsb <- integer(n)
  n1 <- 1:n
  ## Log-likelihood function that is going to be optimized at each step of the FS
  if(one.lambda) {
    # set initial value for optimization when all parameters have to be equal. #
    param <- max(new.lambda[,1])
    loglik <- function(param, data) {
      z.data <- data
      for(kk in new.lambda[, 2]) {
        if(param == 0)
          z.data[, kk] <- log(data[,kk]) * ((exp(mean(log(data[, kk])))))
        else
          z.data[, kk] <- (data[, kk]^param - 1)/(param * ((exp(mean(log(data[, kk])))))^(param - 1))
      }
      det(var(z.data))
    }
  }
  else {
  #Initial values for many values of lambda. #
    param <- new.lambda[, 1]
    loglik <- function(param, data) {
      z.data <- data
      ii <- 0
      for(kk in new.lambda[, 2]){
        ii <- ii+1
        if(param[ii]==0)
          z.data[,kk] <- log(data[, kk]) * ((exp(mean(log(data[, kk])))))
        else
          z.data[, kk] <- (data[, kk]^param[ii] - 1)/(param[ii] * ((exp(mean(log(data[, kk])))))^(param[ii] - 1))
      }
      det(var(z.data))
    }
  }


  eps <- 1e-03
  a <- boundaries[1,1]
  b <- boundaries[1,2]
  transfPar <- function(par){# Ancillary - Forces parameters to be bounded
    transfPar <- rep(0.0, length(par))
    transfPar <- log((par - (a-eps))/((b+eps) - par)) 
    return (transfPar)
   }

   backTransfPar <- function(transfPar){# Ancillary - Forces parameters to be bounded
    par <- rep(0.0, length(transfPar))
    par <- (((b-eps)* exp(transfPar)) + a + eps)/(1+exp(transfPar))
    return (par)
   }

  
   loglik.repar <- function(transfParam, data){#
      param <- backTransfPar(transfParam)
      loglik(param, data)
   }
  message <- 0
  for(m in M:n) {# Start FS.  Optimize over columns of X. #
    Unit[[m-M+1]] <- list(bsb)
    if(length(param) >1){
      fit <- optim(transfPar(param), loglik.repar, data = X[bsb, ])
      param <- backTransfPar(fit$par)
      if(any(param <= (a+eps)) || any(param >= (b-eps))){
        message <- m
        subfit <- 1
        try <- seq(2, 10)
        counter <- 1
        while(subfit > 0){
          a <- boundaries[1,1] * try[counter]
          b <- boundaries[1,2] * try[counter]
          fit <- optim(transfPar(param), loglik.repar, data = X[bsb, ])
          param <- backTransfPar(fit$par)        
          counter <- counter+1
          if(try[counter] == max(try)){
            param <- rep(NA, ncol(Mle[[1]]))
            warning("Algorithm might not have converged in some steps of the search. Boundaries values are treated as missing. Try changing boundaries")
          }
          if(any(param <= (a+eps)) || any(param >= (b-eps))) subfit <- 1
          else subfit <- 0
        }
      }
    }
    else{
      fit <- optimize(loglik, interval = boundaries, data = X[bsb, ])
      param <- fit$minimum
      if(param <= (a+eps) || param >= (b-eps)){
        message <- m
        subfit <- 1
        try <- seq(2, 10)
        counter <- 1
        while(subfit > 0){
          a <- boundaries[1,1] * try[counter]
          b <- boundaries[1,2] * try[counter]
          fit <- optimize(loglik, interval = c(a,b),data = X[bsb, ])
          param <- fit$minimum 
          counter <- counter+1
          if(try[counter] == max(try)){
            param <- rep(NA, ncol(Mle[[1]]))
            warning("Algorithm might not have converged in some steps of the search. Boundaries values are treated as missing. Try changing boundaries")
          }
          if(param <= (a+eps) || param >=(b-eps)) subfit <- 1
          else subfit <- 0
        }
      }

    }
    new.lambda[, 3] <- rep(a, length(new.lambda[, 3]))
    new.lambda[, 4] <- rep(b, length(new.lambda[, 4]))
    

    Mle[[1]][m-M+1, ] <- param
    center <- apply(X.lambda[bsb, ], 2, mean)
    Xm.lambda <- sweep(X.lambda, 2, center)
    R <- qr.R(qr(Xm.lambda[bsb, ]))
    u <- forwardsolve(t(R), t(Xm.lambda))
    distances <- sqrt((m-1) * colSums(u * u))
    names(distances) <- n1
    sorted.distances <- sort(distances)
    if(m < n) {
      newbsb <- as.integer(names(sorted.distances))[1:(m+1)]
      bsb <- newbsb
    }
    
  }
  if (message > 0) cat("Bounds exceeded at step ", message, ". New bounds are ","(",a,",",b,")",sep="","\n")
  value <- list(call = the.call,
                Unit = Unit,
                n = n,
                p = p,
                m = M,
                data = X,
                data.H0 = X.lambda,
                data.names = x.names,
                Mle = Mle,
                H0 = new.lambda,
                forced.onepar = one.lambda,
                dof = length(col.to.transform)
                )
  oldClass(value) <- "fwdtr"
  value
}
