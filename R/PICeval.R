gaussfit = function(pos,width,maxo,PIC){
  x <- PIC[,1]
  y <- PIC[,2]
  a1 <- maxo
  b1 <- pos
  c1 <- 0.25*width
  fit <- nls(y~a1*exp(-((x-b1)/c1)^2),start=list(a1=a1,b1=b1,c1=c1),
             control=list(maxiter=500,tol=1e-3,minFactor=1/512,printEval=F,warnOnly=T))
  yp <- predict(fit)
  sse <- sum((y-yp)^2)
  sst <- sum((y-mean(y))^2)
  r.square <- 1-sse/sst
  return(r.square)
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
  return(sharpness)
}

PICeval = function(PICs){
  options(warn =-1)
  peaks <- PICs$Info
  PICs <- PICs$PICs
  gaussfitness <- c()
  sharpness <- c()
  for (i in 1:nrow(peaks)){
    pos <- peaks[i,'rt']
    width <- max(1,peaks[i,'rtmax']-peaks[i,'rtmin'])
    maxo <- peaks[i,'maxo']
    PIC <- PICs[[i]]
    gaussfitness.i <- 0
    gaussfitness.i <- round(try(gaussfit(pos,width,maxo,PIC)),2)
    sharpness.i <- round(getSharpness(PIC),2)
    gaussfitness <- c(gaussfitness,gaussfitness.i)
    sharpness <- c(sharpness,sharpness.i)
  }
  peaks <- cbind(peaks,gaussfitness,sharpness)
  return(peaks)
}