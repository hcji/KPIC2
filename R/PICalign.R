WhittakerSmooth <- function(y,lambda){
  M <- length(y)
  E <- sparseMatrix(i=1:M,j=1:M,x=1)
  D <- Matrix::diff(E)
  C <- chol(E+lambda*Matrix::t(D)%*%D)
  z <- solve(C,solve(t(C),y))
  return(as.numeric(z))
}

getSpectrum <- function(xset,mz1,mz2,lambda){
  library(Matrix)
  rt_start <- 0
  rt_end <- 10^6
  segSize <- 0
  spectrums <- list()
  s.spectrums <- c()
  for (i in 1:length(xset@path)){
    rt.i <- xset@rt[[i]]
    rt_start <- max(rt_start,rt.i[1])
    rt_end <- min(rt_end,max(rt.i))
    spectrum.i <- rep(0,length(rt.i))
    peak.i <- xset@PICset[[i]]$Info
    id <- which(peak.i[,'mz']>=mz1&peak.i[,'mz']<mz2)
    if (length(id)==0){
      spectrums <- c(spectrums,list(spectrum.i))
      next}
    for (j in 1:length(id)){
      PIC.j <- xset@PICset[[i]]$PICs[[id[j]]]
      spectrum.i[rt.i%in%PIC.j[,1]] <- spectrum.i[rt.i%in%PIC.j[,1]]+PIC.j[,2]
      segSize <- max(segSize,nrow(PIC.j))
    }
    spectrums <- c(spectrums,list(spectrum.i))
  }
  inl <- mean(diff(rt.i))
  rtlist <- seq(rt_start,rt_end,inl)
  for (i in 1:length(xset@path)){
    spectrum.i <- approx(xset@rt[[i]],spectrums[[i]],rtlist)$y
    spectrum.i <- WhittakerSmooth(spectrum.i,lambda)
    s.spectrums <- rbind(s.spectrums,spectrum.i)
  }
  segSize <- max(20,round(2*segSize))
  return(list(rtlist=rtlist,spectra=s.spectrums,segSize=segSize))
}

seg_rtcor <- function(xset,ref,mz1,mz2,shift,lambda){
  r.spectra <- getSpectrum(xset,mz1,mz2,lambda)
  segSize <- r.spectra$segSize
  spectra <- r.spectra$spectra
  if (sum(spectra)<10*length(xset@path)){return(xset)}
  peakmat <- xset@peakmat
  rt.list <- r.spectra$rtlist
  lags <- PAFFT(spectra,spectra[ref,],segSize,shift)$lags
  if (sum(lags)==0){return(xset)}
  for (i in 1:nrow(spectra)){
    s <- findInterval(c(i-1,i),peakmat[,'sample'])
    peak.i <- peakmat[(s[1]+1):s[2],]
    peak.id <- which(peak.i[,'mz']>=mz1&peak.i[,'mz']<=mz2)
    if (length(peak.id)==0){next}
    ids <- (1:length(rt.list))+lags[i,]
    ids[ids<1] <- 1
    ids[ids>length(rt.list)] <- length(rt.list)
    rt.i <- rt.list[ids]
    rb <- max(rt.i[length(rt.i)],rt.list[length(rt.list)])
    rt.list1 <- c(0,rt.list,rb)
    rt.i  <- c(0,rt.i,rb)
    temp <- approx(rt.list1,rt.i,peak.i[peak.id,'rt_cor'])$y
    for (j in 1:length(peak.id)){
      PIC.j <- xset@PICset[[i]]$PICs[[peak.id[j]]]
      rt.dif <- temp[j]-peak.i[peak.id[j],'rt_cor']
      if (length(rt.dif)!=1){stop()}
      PIC.j[,1] <- PIC.j[,1]+rt.dif
      xset@PICset[[i]]$PICs[[peak.id[j]]] <- PIC.j
    }
    peak.i[peak.id,'rt_cor'] <- as.matrix(temp)
    peakmat[(s[1]+1):s[2],] <- peak.i
  }
  xset@peakmat <- peakmat
  return(xset)
}


PAFFT <- function(spectra,reference,segSize,shift){
  lags <- alignedSpectrum <- matrix(0,nrow(spectra),ncol(spectra))
  for (i in 1:nrow(spectra)){
    startpos <- 1
    aligned <- c()
    shifts <- rep(0,ncol(spectra))
    while (startpos <= ncol(spectra)){
      endpos <- startpos+(segSize*2)
      if (endpos>=ncol(spectra)){
        samseg <- spectra[i,startpos:ncol(spectra)]
        refseg <- reference[startpos:ncol(spectra)]
      }else{
        samseg <- spectra[i,(startpos+segSize):(endpos-1)]
        refseg <- reference[(startpos+segSize):(endpos-1)]
        minpos <- findMin(samseg,refseg)
        endpos <- as.numeric(startpos+minpos+segSize)
        samseg <- spectra[i,startpos:endpos]
        refseg <- reference[startpos:endpos]
      }
      lag <- FFTcorr(samseg,refseg,shift)
      samseg <- as.numeric(samseg)
      shifts[startpos:(startpos+length(refseg)-1)] <- lag
      aligned <- c(aligned,as.numeric(move(samseg,lag)))
      startpos <- endpos+1
    }
    alignedSpectrum[i,] <- aligned
    lags[i,] <- shifts
  }
  return(list(alignedSpectrum=alignedSpectrum,lags=lags))
}

FFTcorr <- function(spectrum, target, shift){
  spectrum <- t(as.matrix(spectrum))
  target <- t(as.matrix(target))
  M <- ncol(target)
  diff <- 1000000
  for (i in 1:20){
    curdiff <- ((2^i)-M)
    if (curdiff>0&curdiff<diff){diff <- curdiff}
  }
  target <- cbind(target,t(rep(0,diff)))
  spectrum <- cbind(spectrum,t(rep(0,diff)))
  M <- M+diff
  X <- fft(as.numeric(target))
  Y <- fft(as.numeric(spectrum))
  R <- X*Conj(Y)
  R <- R/M
  rev <- fft(R,inverse=T)/length(rev)
  vals <- Re(rev)
  maxpos <- 1
  maxi <- -1
  if (M<shift){shift <- M}
  for (i in 1:shift){
    if (vals[i] > maxi){
      maxi = vals[i]
      maxpos = i
    }
    if (vals[length(vals)-i+1] > maxi){
      maxi = vals[length(vals)-i+1];
      maxpos = length(vals)-i+1;
    }
  }
  if (maxi < 0.1){lag <- 0
  return(lag)}
  if (maxpos > length(vals)/2){
    lag = maxpos-length(vals)-1
    return(lag)
  }else{lag <- maxpos-1
  return(lag)}
}

move <- function(seg, lag){
  if (lag == 0 || lag >= length(seg)){
    movedSeg <- seg
    return(movedSeg)
  }
  if (lag > 0){
    ins <- rep(1,lag)*seg[1]
    movedSeg <- c(ins,seg[1:(length(seg)-lag)])
    return(movedSeg)
  }
  if (lag < 0){
    lag <- abs(lag);
    ins <- rep(1,lag)*seg[length(seg)]
    movedSeg <- c(seg[(lag+1):length(seg)],ins)
    return(movedSeg)
  }
}

findMin <- function(samseg,refseg){
  Is <- order(samseg)
  Ir <- order(refseg)
  minposA <- c()
  minInt <- c()
  minpos <- NA
  for (i in 1:round(length(Is)/20)){
    for (j in 1:round(length(Is)/20)){
      if (Ir[j]==Is[i]){minpos <- Is[i]
      return(minpos)}
    }
  }
  if (is.na(minpos)){minpos <- Is[1]}
  return(minpos)
}