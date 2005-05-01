fwdtrProfilePlot <- function(x, psfrag.labels = FALSE)
{
  n <- x$n
  m <- x$step.fwd
  prof <- x$profile
  number.plot <- x$p
  xaxis <- x$lambda
  ci <- x$ci
  onepar <- x$onepar

	if(!interactive())
	  return()

  temp.range <- lapply(prof, range)

  temp.low1 <- max(unlist(lapply(temp.range, min))) 
  temp.low2 <- max(unlist(lapply(temp.range, max))) -20
  temp.low <- min(temp.low1, temp.low2)
  temp.high <- max(unlist(lapply(temp.range, max)))
  
  ylim <- c(temp.low, temp.high)
  if(onepar){
    if(psfrag.labels) {
      x.label <- "xlab"
      y.label <- "ylab"
      main.label <- "main"
    }
    else{
      x.label <- "lambda"
      y.label <- "log-likelihood"
      main.label <- " "
    }
    plot(xaxis,
         prof[[1]],
         type="l",
         xlab = x.label,
         ylab = y.label,
         main = main.label,
         ylim = ylim)
    abline(v = ci[[1]], col = 8)
    points(x$Mle, ylim[1], pch=20)
    text(xaxis[round(length(xaxis)/6)], ylim[2], paste("m = ", m, sep=""))
    
  }
  
  
  else{
    choices <- x$x.names
    tmenu <- paste("plot:", choices)
    if(length(choices) > 1) tmenu <- c(tmenu, "plot: all")
    pick <- 1
    while(pick > 0) {
      pick <- menu(tmenu, title =
                   "\nMake a plot selection (or 0 to exit):")
      
      if(pick >0){
        if(pick < length(tmenu)){
          if(psfrag.labels) {
            x.label <- paste("xlab",pick,sep="")
            y.label <- paste("ylab",pick,sep="")
            main.label <- paste("main",pick,sep="")
          }
          else{
            x.label <- "lambda"
            y.label <- x$x.names[pick]
            main.label <- " "
          }
          
          plot(xaxis,
               prof[[pick]],
               type="l",
               xlab = x.label,
               ylab = y.label,
               main = main.label,
               ylim = ylim)
          abline(v = ci[[pick]], col = 8)
          points(x$Mle[pick], ylim[1], pch=20)
          text(xaxis[round(length(xaxis)/6)], ylim[2], paste("m = ", m, sep=""))
        }
        if(pick == length(tmenu)) {
          oldpar <- par(no.readonly = TRUE)  
          if (number.plot <=8){
            a <- 2
            b <- round(number.plot/2)+1
            if((b-a) > 1) b <- b-1
          } 
          else {
            a <- 3
            b <- round(number.plot/3)+1
          }
          if(number.plot == 1) a <- b <- 1
          if(number.plot == 2) b <- 1
          par(mfrow = c(a,b))
          for (i in 1:number.plot){
            if(psfrag.labels) {
              x.label <- paste("xlab",i,sep="")
              y.label <- paste("ylab",i,sep="")
              main.label <- paste("main",i,sep="")
            }
            else{
              x.label <- "lambda"
              y.label <- x$x.names[i]
              main.label <- " "
            }
            plot(xaxis, prof[[i]], type="l", col = 1, xlab = x.label, ylab = y.label, main = main.label, ylim = ylim)
            abline(v = ci[[i]], col = 8)
            points(x$Mle[i], ylim[1], pch=20)
            text(xaxis[round(length(xaxis)/6)], ylim[2], paste("m = ", m, sep=""))
          }
          par(oldpar)
        }
      }
    }
  }
  invisible()
}


