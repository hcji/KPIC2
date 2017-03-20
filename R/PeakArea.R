cwt_haar <- function(ms, scales=1) {
  psi_xval <- seq(0,1,length=1024)
  psi <- c(0,rep(1,511),rep(-1,511),0)

  oldLen <- length(ms)
  ms <- extendNBase(ms, nLevel=NULL, base=2)
  len <- length(ms)
  nbscales <- length(scales)
  wCoefs <- NULL
  
  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax  <- psi_xval[length(psi_xval)]
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    if (length(j) == 1)		j <- c(1, 1)
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
    wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
    ## Shift the position with half wavelet width
    wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
    wCoefs <- cbind(wCoefs, wCoefs.i)
  }
  if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
  colnames(wCoefs) <- scales
  wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
  return(wCoefs)
}

localMaximum <- function (x, winSize = 5) {
    len <- length(x)
    rNum <- ceiling(len/winSize)
    
    ## Transform the vector as a matrix with column length equals winSize
    ##		and find the maximum position at each row.
    y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    
    ## keep the result
    localMax <- rep(0, len)
    localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1
    
    ## Shift the vector with winSize/2 and do the same operation
    shift <- floor(winSize/2)
    rNum <- ceiling((len + shift)/winSize)	
    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1
    
    ## Check whether there is some local maxima have in between distance less than winSize
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
      selMaxInd1 <- maxInd[selInd]
      selMaxInd2 <- maxInd[selInd + 1]
      temp <- x[selMaxInd1] - x[selMaxInd2]
      localMax[selMaxInd1[temp <= 0]] <- 0
      localMax[selMaxInd2[temp > 0]] <- 0
    }
    
    return(localMax)
}

extendNBase <- function(x, nLevel=1, base=2, ...) {
  if (!is.matrix(x)) {
    x <- matrix(x, ncol=1)
  } else if (min(dim(x)) == 1) {
    x <- matrix(x, ncol=1)
  }
  
  nR <- nrow(x)
  if (is.null(nLevel)) {
    nR1 <- nextn(nR, base)		
  } else {
    nR1 <- ceiling(nR / base^nLevel) * base^nLevel		
  }
  if (nR != nR1) {
    x <- extendLength(x, addLength=nR1-nR, ...)
  }
  
  return(x)
}

extendLength <- function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) {
  if (is.null(addLength)) stop('Please provide the length to be added!')
  if (!is.matrix(x)) x <- matrix(x, ncol=1)	
  method <- match.arg(method)
  direction <- match.arg(direction)
  
  nR <- nrow(x)
  nR1 <- nR + addLength
  if (direction == 'both') {
    left <- right <- addLength
  } else if (direction == 'right') {
    left <- 0
    right <- addLength
  } else if (direction == 'left') {
    left <- addLength
    right <- 0
  }
  
  if (right > 0) {
    x <- switch(method,
                reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
                open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
                circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
  }
  
  if (left > 0) {
    x <- switch(method,
                reflection =rbind(x[addLength:1, , drop=FALSE], x),
                open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
                circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
  }
  if (ncol(x) == 1)  x <- as.vector(x)
  
  return(x)
}

widthEstimationCWT <- function(x,majorPeakInfo) {
  wCoefs_haar <- cwt_haar(x, 1:ceiling(max(majorPeakInfo$peakScale)))
  peakIndex <- majorPeakInfo$peakIndex
  peakScale <- majorPeakInfo$peakScale
  
  peakWidth <- list()
  peakWidth[["peakIndex"]]= majorPeakInfo$peakIndex
  peakWidth[["peakIndexLower"]]= majorPeakInfo$peakIndex
  peakWidth[["peakIndexUpper"]]= majorPeakInfo$peakIndex
  
  for(i in 1:length(peakIndex)){
    peakIndex.i=peakIndex[i]
    peakScale.i=round(peakScale[i])
    wCoefs_haar.i=wCoefs_haar[,peakScale.i]
    wCoefs_haar.i.abs=abs(wCoefs_haar.i)
    localmax=localMaximum(-wCoefs_haar.i.abs,winSize=5)      
    #localmax=localmax & abs(wCoefs_haar.i)<(mean(wCoefs_haar.i.abs[localmax==1])+0.5*sd(wCoefs_haar.i.abs[localmax==1]))
    localmax=as.numeric(localmax)
    localmax[peakIndex]=0
    temp1 <- max(1,(peakIndex.i-peakScale.i/2+1))
    temp2 <- min(length(localmax),(peakIndex.i+peakScale.i/2-1))
    localmax[temp1:temp2]=0
    
    Lef =0
    Rig =0
    
    peakScale.i.3=2*peakScale.i
    
    if(i==1){
      maxIndexL=max(c((peakIndex.i-peakScale.i.3),1))
    }else{
      maxIndexL=max(c((peakIndex.i-peakScale.i.3),peakIndex[i-1]))
    }
    
    if(i==length(peakIndex)){
      minIndexR=min(c((peakIndex.i+peakScale.i.3),length(localmax)))
    } else{
      minIndexR=min(c((peakIndex.i+peakScale.i.3),peakIndex[i+1]))
    }
    ignoreL=1:maxIndexL
    ignoreR=minIndexR:length(localmax)        
    localmax[c(ignoreL,ignoreR)]=0
    
    localmax[c(peakIndex.i,temp1:temp2)]=0
    
    bi=which(localmax==1)
    
    biLeft=bi[bi<peakIndex.i]
    useL= maxIndexL:peakIndex.i
    minIndexLeft=useL[which(min(x[useL])==x[useL])]
    
    if(length(biLeft)==0){
      Lef=minIndexLeft
    }else{
      minbaselineIndexLeft=biLeft[which(min(x[biLeft])==x[biLeft])]
      if(minIndexLeft>=(peakIndex.i-peakScale.i/2+1)){
        Lef=minbaselineIndexLeft
      }else{
        Lef=max(c(minIndexLeft,minbaselineIndexLeft))
      }   
    }
    
    biRight=bi[bi>peakIndex.i]
    useR= peakIndex.i:minIndexR
    minIndexRight=useR[which(min(x[useR])==x[useR])]
    
    if(length(biRight)==0){
      Rig= minIndexRight
    }else{
      minbaselineIndexRight=biRight[which(min(x[biRight])==x[biRight])] 
      if(minIndexRight<=(peakIndex.i+peakScale.i/2-1)){
        Rig=minbaselineIndexRight
      }else{
        Rig=min(c(minIndexRight,minbaselineIndexRight))
      }      
    }
    # peakWidth[[paste(peakIndex.i)]]=Lef:Rig
    peakWidth[["peakIndexLower"]][i]= Lef
    peakWidth[["peakIndexUpper"]][i]= Rig
  }      
  return(peakWidth)		
}

integration <- function(x,yf){
  # The yf is a collection of each fitting peaks by calling the area function.
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }
  
  # check if the integrand is numeric
  if (!is.numeric(yf))
  {
    stop('The integrand "yf" is not numeric.')
  }
  
  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(yf))
  {
    stop('The lengths of the variable of integration and the integrand do not match.')
  }
  
  ### finish checks
  
  # obtain length of variable of integration and integrand
  n = length(x)
  
  # integrate using the trapezoidal rule
  integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  
  # print the definite integral
  return(integral)
}

