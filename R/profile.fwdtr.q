profile.fwdtr <- function(fitted, step.fwd = NULL, conf = 0.95, bounds = NULL, ...)
{
  n <- fitted$n
  H0 <- fitted$H0
  p <- fitted$p
  col.maximized <- H0[, 2]
  m <- fitted$m
  if(is.null(bounds)) bounds <- fitted$H0[1, 3:4]
  if((bounds[1] !=fitted$H0[1, 3]) || (bounds[2] !=fitted$H0[1, 4]))
      cat("Profile likelihood might show maxima different from constrained optimazion", "\n")
  Mle <- matrix(NA, nrow(fitted$Mle[[1]]), p)
  onepar <- fitted$forced.onepar
  jj <- 0
  if(onepar == TRUE) Mle <- fitted$Mle[[1]]
  else{
    for(i in col.maximized){
      jj <- jj+1
      Mle[, i] <- fitted$Mle[[1]][, jj]
    }
  }
  x.names <- fitted$data.names[[2]][col.maximized]
  if(is.null(step.fwd)) step.fwd <- n
  if (step.fwd > n) stop("Step of the forward search cannot exceed the number of units")
  if(step.fwd < m) stop(paste("The Search starts at step", m,sep=" "))
  the.call <- match.call()
  if(is.null(x.names))
    x.names <- list(1:length(col.maxmized), paste("V",col.maximized, sep = ""))

  bsb <- fitted$Unit[[step.fwd-m+1]][[1]] # Get data from the specified step of the search
  max.mle <- Mle[(step.fwd-m+1), ]
  temp.data <- as.matrix(fitted$data)
  
  my.bsb <- temp.data[bsb, ]
  prof <- list()
  ci <- list()
  profloglik <- function(param, col.param){
    z.data <- my.bsb
    mmm <- max(my.bsb)
    myscale <- apply(my.bsb, 2, max)
    my.bsb <- sweep(my.bsb, 2, myscale,FUN = "/")
    if(param == 0){
      z.data[, col.param] <- log(my.bsb[, col.param]) * ((exp(mean(log(my.bsb[, col.param])))))
      }
    else{
      z.data[, col.param] <- (my.bsb[, col.param]^param - 1)/(param * ((exp(mean(log(my.bsb[, col.param])))))^(param - 1))
    	}
    jj <- 0
    for(kk in col.fixed) {
      jj <- jj+1
      if(par.fixed[jj]==0){
        z.data[, kk] <- log(my.bsb[, kk]) * ((exp(mean(log(my.bsb[, kk])))))
        }
      else{
        z.data[, kk] <- (my.bsb[, kk]^par.fixed[jj] - 1)/(par.fixed[jj] * ((exp(mean(log(my.bsb[, kk])))))^(par.fixed[jj] - 1))
      }
      }
    det(var(z.data*mmm))
  }
  
  if(onepar == TRUE) {
     loglik <- function(param, col.maximized) {
      z.data <-  my.bsb
      for(kk in col.maximized) {
        if(param == 0)
          z.data[, kk] <- log(my.bsb[,kk]) * ((exp(mean(log(my.bsb[, kk])))))
        else
          z.data[, kk] <- (my.bsb[, kk]^param - 1)/(param * ((exp(mean(log(my.bsb[, kk])))))^(param - 1))
      }
      det(var(z.data))
    }
     lambda <- seq(bounds[1], bounds[2], 0.1)
     prof[[1]] <- rep(NA, length(lambda))
     i <- 0
     j <- 0

     for(j in 1:length(lambda)){
       i <- i+1

       prof[[1]][j] <- -step.fwd*log(loglik(lambda[i], col.maximized = col.maximized))
     }

     xtemp <- lambda[prof[[1]] < (max(prof[[1]])-qchisq(conf, df=1)/2)]
     if(length(max(diff(xtemp))) > 1) stop("Increase the width of your bands")

     index <- which(diff(xtemp) == max(diff(xtemp)))
     index <- c(index, index+1)
     ci[[1]] <- xtemp[index]
     
   }
  else{
    message <- 0
      jjj <- 0
      for(ll in col.maximized) {
        jjj <- jjj+1
        col.fixed <- col.maximized[-jjj]
        par.fixed <- Mle[(step.fwd-m+1), col.fixed]
        par.fixed[is.na(par.fixed)] <- 1
        lambda <- seq(bounds[1], bounds[2], 0.1)
        prof[[jjj]] <- rep(NA, length(lambda))
        i <- 0
        for(j in 1:length(lambda)){
          i <- i+1
          prof[[jjj]][j] <- -step.fwd*log(profloglik(lambda[i], col.param = ll))
        }
        xtemp <- lambda[prof[[jjj]] < (max(prof[[jjj]])-qchisq(conf, df=1)/2)]
        cat(length(xtemp),"\n")
        if(length(xtemp)<1) {
            ci[[jjj]] <- NULL
            message <- message +1
        }
        else {
        if((max(diff(xtemp)) - max(diff(lambda))) <=0){
          if(((min(xtemp) < 0) && (max(xtemp) <0)) && (max.mle[jjj] > 0)) ci[[jjj]] <- max(xtemp)
          if(((min(xtemp) < 0) && (max(xtemp) <0)) && (max.mle[jjj] < 0)) ci[[jjj]] <- min(xtemp)
          if(((min(xtemp) < 0) && (max(xtemp) > 0)) && (max.mle[jjj] <0)) ci[[jjj]] <- min(xtemp)
          if(((min(xtemp) < 0) && (max(xtemp) > 0)) && (max.mle[jjj] > 0)) ci[[jjj]] <- max(xtemp)
          if((min(xtemp) > 0) && (max(xtemp) > 0)) ci[[jjj]] <- min(xtemp)
          message <- message +1
        }
        else { 
          index <- which(diff(xtemp) == max(diff(xtemp)))
          index <- c(index, index+1)
          ci[[jjj]] <- xtemp[index]
        }
      }
    }
    if (message > 0) cat("Confidence bands approximated: flat likelihood or MLE close to bounds", "\n")
  }
  if(onepar){
    names(ci) <- "common.parameter"
    names(max.mle) <- "common.parameter"
  }
  else{
    names(ci) <- x.names
    names(max.mle[!is.na(max.mle)]) <- x.names
  }
  out <- list(call = the.call,
              lambda=lambda,
              profile=prof,
              ci = ci,
              x.names = x.names,
              step.fwd = step.fwd,
              onepar = onepar,
              p = length(col.maximized), 
              Mle = max.mle[!is.na(max.mle)]
              )
  oldClass(out) <- "profile.fwdtr"
  out
}



