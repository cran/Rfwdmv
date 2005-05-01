plot.fwdtr.test <- function(x, psfrag.labels = FALSE, ...)
{
  n <- x$n
  m <- x$m
  test <- x$test
  lambda.around <- x$lambda.around
  n.lambda <- length(lambda.around)
  number.plot <- length(x$col.to.compare)

  if(!interactive())
	  return()

  oldpar <- par(no.readonly = TRUE)
  choices <- x$col.names[x$col.to.compare]
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
            x.label <- "Subset Size"
            y.label <- x$col.names[x$col.to.compare[pick]]
            main.label <- " "
          }
          
          matplot(m:n, test[[pick]], type="l", col = 1, lty = 1:n.lambda, lwd = 1:n.lambda, xlab = x.label, ylab = y.label, main = main.label)
          abline(h = -2.58, col = 8)
          abline(h = 2.58, col = 8)
          labels <- as.character(lambda.around)
          if(any(abs(lambda.around - 1/3) < 0.0000000001)) {
            labels[abs(lambda.around - 1/3) < 0.0000000001] <- "1/3"
          }
          labels <- paste(" ", labels, sep = "")
          text(n, test[[pick]][nrow(test[[pick]]), ], labels, adj = 0)
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
              x.label <- "Subset Size"
              y.label <- x$col.names[x$col.to.compare[i]]
              main.label <- " "
            }
            matplot(m:n, test[[i]], type="l", col = 1, lty = 1:n.lambda, lwd = 1:n.lambda, xlab = x.label, ylab = y.label, main = main.label)
            abline(h = -2.58, col = 8)
            abline(h = 2.58, col = 8)
            labels <- as.character(lambda.around)
            if(any(abs(lambda.around - 1/3) < 0.0000000001)) {
              labels[abs(lambda.around - 1/3) < 0.0000000001] <- "1/3"
            }
            labels <- paste(" ", labels, sep = "")
            text(n, test[[i]][nrow(test[[i]]), ], labels, adj = 0)
          }
        }
        par(oldpar)
      }
  }
  
  invisible()
}
