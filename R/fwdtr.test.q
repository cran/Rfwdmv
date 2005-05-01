fwdtr.test <- function(X, parameters, n.bsb = 50, col.to.compare = "all", lambda.around = c(-1.0, -0.5, 0.0, 0.5, 1.0),
                              one.lambda = FALSE)
{
  X <- as.matrix(X)
  X.names <- dimnames(X)
  the.call <- match.call()
  X.parameters <- X
  if(is.character(col.to.compare)) {
    col.to.compare <- seq(1:ncol(X))
  }
  
  else {
    col.to.compare <- as.integer(col.to.compare)
  }
  
  p <- dim(X)[2]
  number.of.na <- 0
  message <- NULL
  
  compute.likelihood.ratio.test <- function(x){
    p <- x$p
    n <- x$n
    m <- x$m
    dof <- x$dof
    n.steps <- n-m+1
    constraint <- rep(1.0, p)
    constraint[x$H0[, 2]] <- x$H0[, 1]
    if(x$forced.onepar) {
      constraint <- rep(min(x$H0[,1]), p)
      x$Mle[[1]] <- matrix(x$Mle[[1]], ncol=p, nrow = nrow(x$Mle[[1]]))
    }
    likratio <- NULL
    data <- x$data
    for(i in 1:n.steps) {
      sub.data <- data[x$Unit[[i]][[1]], ]
      for (ii in 1:length(constraint)) {
        G <- exp(mean(log(sub.data[,ii])))
        if(constraint[ii] == 0) {
          sub.data[, ii] <- log(sub.data[, ii])*G
        }
        else {
          sub.data[,ii] <- ((sub.data[,ii]^constraint[ii]) - 1)/(constraint[ii] * G^(constraint[ii]-1))
        }
      }
      sub.data.H0 <- sub.data
      sub.data[,x$H0[,2]] <- data[x$Unit[[i]][[1]], x$H0[, 2]]
      j <- 0
      for(k in x$H0[,2]) {
        j <- j+1
        sub.data[, k] <- (sub.data[,k]^x$Mle[[1]][i,j] - 1)/(x$Mle[[1]][i,j] * ((exp(mean(log(sub.data[,k])))))^(x$Mle[[1]][i,j] - 1))
      }
      likratio[i] <- nrow(sub.data) * log(det(var(sub.data.H0))/det(var(sub.data)))
      if(likratio[i] < 0){
        likratio[i] <- NA
      }
    }
    likratio
  }

  if(is.null(X.names))
    X.names <- list(1:dim(X)[1], paste("V", 1:p, sep = ""))
  
  if(length(parameters) != p)
    stop("the argument \"parameters\" must be a vector of p elements.")
  lik.ratio.test <- list()
  
  n.bsb <- max((p+1), as.integer((n.bsb * nrow(X))/100))
  kk <- 0
  for (j in col.to.compare) {
 #   cat("Expanding each lambda around = ",parameters[j],sep="","\n")    
    kk <- kk + 1
    X.parameters <- X
    for (ii in 1:length(parameters)){
      G <- exp(mean(log(X.parameters[,ii])))
      if(parameters[ii] == 0) X.parameters[, ii] <- log(X[, ii])*G
      else{
        X.parameters[,ii] <- ((X[,ii]^parameters[ii]) - 1)/(parameters[ii] * G^(parameters[ii]-1))
      }
    }
    
    Xj <- X.parameters
    Xj[,col.to.compare[kk]] <- X[, col.to.compare[kk]]
    for(i in 1:length(lambda.around)){
      cat("Expanding lambda = ",lambda.around[i], " around "  ,parameters[j], ". ",sep ="","\n")
      if(one.lambda) temp <- fwdtr(Xj, n.bsb = n.bsb, lambda =lambda.around[i], col.to.transform = col.to.compare[kk], one.lambda = TRUE)
      else temp <- fwdtr(Xj, n.bsb = n.bsb, lambda =lambda.around[i], col.to.transform = col.to.compare[kk])
      
      n <- temp$n
      m <- temp$m
      if(i ==1) lik.ratio.test[[kk]] <- matrix(as.double(NA), nrow =  n-m+1, ncol = length(lambda.around))
      lik.ratio.test[[kk]][1:(n-m+1), i] <- sign(temp$Mle[[1]] - lambda.around[i])*sqrt(compute.likelihood.ratio.test(temp))
      number.of.na <- sum(c(number.of.na, sum(is.na(compute.likelihood.ratio.test(temp)))))
    }
  }
  if(number.of.na >0) message <- paste("optim failed to converge ", number.of.na, " times",sep="")
  out <- list(call = the.call,
              test = lik.ratio.test,
              n=n,
              m=m,
              col.names = X.names[[2]],
              col.to.compare = col.to.compare,
              lambda.around = lambda.around,
              message = message)
  oldClass(out) <- "fwdtr.test"
  out
}
