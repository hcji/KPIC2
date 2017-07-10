gaussfit = function(PIC){
  x <- PIC[,1]
  y <- PIC[,2]
  a1 <- max(y)
  b1 <- x[which.max(y)]
  c1 <- 0.25*length(x)
  res <- 0
  res <- try({
    fit <- nls(y~a1*exp(-((x-b1)/c1)^2),start=list(a1=a1,b1=b1,c1=c1),
               control=list(maxiter=5000,tol=1e-2,minFactor=1/512,printEval=F,warnOnly=T))
    yp <- predict(fit)
    sse <- sum((y-yp)^2)
    sst <- sum((y-mean(y))^2)
    r.square <- 1-sse/sst
    res <- r.square
  })
  return(round(res,2))
}

getSharpness = function(PIC){
  sig <- PIC[,2]
  sig <- c(sig[1]-0.001,sig,sig[length(sig)]-0.001)
  sig <- sig+abs(rnorm(length(sig)))/10000
  p <- which(sig==max(sig))
  sharpness <- 0
  for (i in 2:p){
    temp <- ((sig[i]-sig[i-1]))/sig[i-1]
    sharpness <- sharpness+temp
  }
  for (j in p:length(sig)){
    temp <- ((sig[i]-sig[i+1])/sig[i+1])
    sharpness <- sharpness+temp
  }
  return(round(sharpness,2))
}
