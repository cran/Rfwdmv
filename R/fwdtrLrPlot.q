fwdtrLrPlot <- function(x, psfrag.labels = FALSE)
{
  p <- x$p
  n <- x$n
  m <- x$m
  dof <- x$dof
  n.steps <- n-m+1
  onepar <- x$forced.onepar
  
  if(psfrag.labels) {
    x.label <- "xlab"
    y.label <- "ylab"
    main.label <- "main"
  }
  else{
    x.label <- "Subset Size"
    y.label <- "Likelihood Ratio"
    main.label <- " "
  }


  
  constraint <- rep(1.0, p)
  
  constraint[x$H0[, 2]] <- x$H0[, 1]
  
  if(onepar) {  # Loglik ratio when all paramters are forced to be equal in the numerical maximization step. #
    constraint <- rep(min(x$H0[,1]), p)
    x$Mle[[1]] <- matrix(x$Mle[[1]], ncol=p, nrow = nrow(x$Mle[[1]]))
    dof <- 1 # Chi square 1 for confidence bands
  }
  L.max.Mle <- NULL

  # Read the data from original matrix. #
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
    
    L.max.Mle[i] <- nrow(sub.data) * log(det(var(sub.data.H0))/det(var(sub.data)))
  }
  
  plot(m:n, L.max.Mle, type = "l", xlab = "Subset Size", ylab = "Likelihood Ratio")
  
  abline(h = qchisq(0.95, dof), col = 8)
  abline(h = qchisq(0.99, dof), col = 8)

  invisible()
}

