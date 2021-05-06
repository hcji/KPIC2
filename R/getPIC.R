getPIC <- function(raw, level, mztol=0.1, gap=3, width=5, min_snr=4, export='FALSE',...){
  orders <- order(raw$mzs)
  mzs <- raw$mzs[orders]
  scans <- raw$scans[orders]
  ints <- raw$ints[orders]
  scantime <- raw$times
  path <- raw$path
  rm(raw)

  # set seeds
  scanwidth <- as.integer(width/mean(diff(scantime[scantime>0])))
  seeds <- which(ints>level)
  seeds <- seeds[order(-ints[seeds])]
  clu <- rep(0, length(ints))

  # detect pics
  clu <- getPIP(seeds,scans,mzs,clu,mztol,gap)
  orders <- order(clu)
  clu <- clu[orders]
  mzs <- mzs[orders]
  scans <- scans[orders]
  ints <- ints[orders]

  picind <- c(findInterval(0:max(clu),clu))
  pics <- lapply(1:(length(picind)-1),function(s){
    pic <- cbind(scans[(picind[s]+1):picind[s+1]],ints[(picind[s]+1):picind[s+1]],mzs[(picind[s]+1):picind[s+1]])
    pic <- pic[order(pic[,1]),]
    return(pic)
  })
  pic_length <- unlist(sapply(pics,length))
  pics <- pics[pic_length>3*scanwidth]

  # interpolation of missing points of pic
  pics <- lapply(pics,function(pic){
    scan <- pic[1,1]:pic[nrow(pic),1]
    int <- approx(pic[,1],pic[,2],scan)$y
    mz <- approx(pic[,1],pic[,3],scan)$y
    cbind(scan,int,mz)
  })
  gc()

  # peak detection
  peaks <- lapply(pics,function(pic){
    peak_detection(pic[,2], min_snr, level)
  })
  nps <- sapply(peaks,function(peaki){
    length(peaki$peakIndex)
  })
  pics <- pics[nps>0]
  peaks <- peaks[nps>0]
  gc()

  output <- list(path=path, scantime=scantime, pics=pics, peaks=peaks)
  if (export){
    exportJSON <- toJSON(output)
    splitname <- strsplit(path,"\\.")[[1]][1]
    outpath <- paste(splitname,'json',sep='.')
    write(exportJSON,outpath)
    gc()
  }
  return(output)
}

getPIC.kmeans <- function(raw, level, mztol=0.1, gap=3, width=c(5,60), alpha=0.3, min_snr=4, export='FALSE', ...){

  orders <- order(raw$mzs)
  mzs <- raw$mzs[orders]
  scans <- raw$scans[orders]
  ints <- raw$ints[orders]
  notused <- rep(TRUE, length(ints))
  scantime <- raw$times
  path <- raw$path
  rm(raw)

  # set seeds
  scanwidth <- as.integer(width/mean(diff(scantime[scantime>0])))
  min_width <- scanwidth[1]
  max_width <- scanwidth[2]
  seeds <- which(ints>level)
  seeds <- seeds[order(-ints[seeds])]

  # detect pics
  pics <- list()
  for (seed in seeds){
    if (!notused[seed]){next}
    ref_mz <- mzs[seed]
    ref_scan <- scans[seed]
    ref_int <- ints[seed]

    roi <- 1+getROI(seed, scans, mzs, ints, notused, mztol, max_width)
    roi_mzs <- mzs[roi]

    diffs <- (roi_mzs-ref_mz)^2
    clu.res <- Ckmeans.1d.dp(diffs, c(1,5))
    sel_id <- roi[which(clu.res$cluster==which.min(clu.res$centers))]
    sel_mz <- mzs[sel_id]
    sel_scan <- scans[sel_id]
    sel_ints <- ints[sel_id]

    pic_id <- unique(collectPIC(ref_scan, ref_mz, ref_int, sel_id, sel_scan, sel_mz, sel_ints, gap, alpha))
    notused[pic_id] <- FALSE
    pic <- data.frame(scans[pic_id], ints[pic_id], mzs[pic_id])
    pics <- c(pics, list(pic))
  }
  gc()

  pic_length <- unlist(sapply(pics,nrow))
  pics <- pics[pic_length>scanwidth[1]]

  # interpolation of missing points of pic
  pics <- lapply(pics,function(pic){
    scan <- min(pic[,1]):max(pic[,1])
    int <- approx(pic[,1],pic[,2],scan)$y
    mz <- pic[match(scan, pic[,1], nomatch = NA),3]
    cbind(scan,int,mz)
  })
  gc()

  # peak detection
  peaks <- lapply(pics,function(pic){
    peak_detection(pic[,2], min_snr, level)
  })
  nps <- sapply(peaks,function(peaki){
    length(peaki$peakIndex)
  })
  pics <- pics[nps>0]
  peaks <- peaks[nps>0]
  gc()

  output <- list(path=path, scantime=scantime, pics=pics, peaks=peaks)
  if (export){
    exportJSON <- toJSON(output)
    splitname <- strsplit(filename,"\\.")[[1]][1]
    outpath <- paste(splitname,'json',sep='.')
    write(exportJSON,outpath)
  }
  return(output)
}

.PICsplit <- function(peak,pic){
  valls <- sapply(1:(length(peak$peakIndex)-1), function(s){
    peak$peakIndex[s]-1+which.min(pic[peak$peakIndex[s]:peak$peakIndex[s+1],2])
  })
  vall.vals <- pic[valls,2]
  thres <- sapply(1:(length(peak$peakIndex)-1), function(s){
    0.5*min(pic[c(peak$peakIndex[s],peak$peakIndex[s+1]),2])
  })
  sp <- which(vall.vals < thres)
  sp <- sp[diff(c(0,valls[sp]))>10]
  if (length(sp)>0){
    splp <- valls[sp]
    starts <- c(1,splp+1)
    ends <- c(splp,nrow(pic))
    pic <- lapply(1:length(starts),function(s){
      pic[starts[s]:ends[s],]
    })
    ps <- c(1,sp+1)
    pe <- c(sp,length(peak$peakIndex))
    peak <- lapply(1:length(ps),function(s){
      peakIndex <- peak$peakIndex[ps[s]:pe[s]]-starts[s]+1
      snr <- peak$snr[ps[s]:pe[s]]
      signals <- peak$signals[ps[s]:pe[s]]
      peakScale <- peak$peakScale[ps[s]:pe[s]]
      list(peakIndex=peakIndex,snr=snr,signals=signals,peakScale=peakScale)
    })
    return(list(n=length(sp),peak=peak,pic=pic))
  } else {return(list(n=0,peak=peak,pic=pic))}
}

PICsplit <- function(pics){
  pics1 <- list()
  peaks1 <- list()
  for (i in 1:length(pics$pics)){
    peak <- pics$peaks[[i]]
    pic <- pics$pics[[i]]
    if (length(pics$peaks[[i]]$peakIndex)>1){
      res <- .PICsplit(peak,pic)
      if (res$n>0){
        pics1 <- c(pics1,res$pic)
        peaks1 <- c(peaks1,res$peak)
        next
      }
    }
    pics1 <- c(pics1,list(pic))
    peaks1 <- c(peaks1,list(peak))
  }

  pics$pics <- pics1
  pics$peaks <- peaks1
  return(pics)
}

getPeaks <- function(pics){
  mzinfo <- lapply(pics$pics,function(pic){
    mz <- mean(pic[,3], na.rm=TRUE)
    mzmin <- min(pic[,3], na.rm=TRUE)
    mzmax <- max(pic[,3], na.rm=TRUE)
    mzrsd <- sd(pic[,3], na.rm=TRUE)/mz*10^6
    c(mz,mzmin,mzmax,mzrsd)
  })

  if (!is.null(pics$peaks)){
    peakpos <- sapply(pics$peaks,function(peaki){
      peaki$peakIndex[which.max(peaki$signals)]
    })
    snr <- sapply(pics$peaks,function(peaki){
      peaki$snr[which.max(peaki$signals)]
    })
    snr <- round(snr,2)
  } else {
    peakpos <- sapply(pics$pics,function(pic){
      pic[which.max(pic[,2]),1]
    })
    snr <- rep(NA, length(peakpos))
  }

  rt <- sapply(1:length(pics$pics),function(s){
    pics$scantime[pics$pics[[s]][peakpos[s],1]]
  })

  maxo <- sapply(pics$pics,function(pic){
    max(pic[,2])
  })

  rtmin <- sapply(pics$pics,function(pic){
    pics$scantime[pic[1,1]]
  })
  rtmax <- sapply(pics$pics,function(pic){
    pics$scantime[pic[nrow(pic),1]]
  })

  area <- sapply(pics$pics,function(pic){
    round(integration(pics$scantime[pic[,1]],pic[,2]))
  })

  mzinfo <- round(do.call(rbind,mzinfo),4)
  colnames(mzinfo) <- c('mz','mzmin','mzmax','mzrsd')

  peakinfo <- cbind(rt,rtmin,rtmax,mzinfo,maxo,area,snr)
  pics$peakinfo <- peakinfo

  return(pics)
}

.PICresolve <- function(peak, pic, pval){
  hh <- 0.5*(pic[c(peak$peakIndex),2])
  ranges <- whichAsIRanges(pic[,2]>min(hh))

  mids <- round(diff(peak$peakIndex)/2) + peak$peakIndex[1:(length(peak$peakIndex)-1)]
  starts <- c(min(start(ranges)),mids)
  ends <- c(mids,max(end(ranges)))

  tpeak <- rep(TRUE,length(starts))

  for (i in 1:(length(starts)-1)){
    mz1 <- pic[starts[i]:ends[i],3]
    mz2 <- pic[starts[i+1]:ends[i+1],3]
    if (min(length(mz1),length(mz2))<2) {
      p <- 1
    } else {
      p <- t.test(mz1,mz2)$p.value
    }
    if (p > pval){
      if (peak$signals[i]>peak$signals[i+1]){tpeak[i+1] <- FALSE
      }else{tpeak[i] <- FALSE}
    }
  }

  peak$peakIndex <- peak$peakIndex[tpeak]
  peak$snr <- peak$snr[tpeak]
  peak$signals <- peak$signals[tpeak]
  peak$peakScale <- peak$peakScale[tpeak]

  return(list(peak=peak, pic=pic))
}

PICresolve <- function(pics, pval=0.01){

  pics1 <- list()
  peaks1 <- list()

  npeak <- sapply(pics$peaks,function(peak){
    length(peak$peakIndex)
  })

  for (i in 1:length(pics$pics)){
    pic <- pics$pics[[i]]
    peak <- pics$peaks[[i]]
    if (npeak[i] < 2){
      pics1 <- c(pics1,list(pic))
      peaks1 <- c(peaks1,list(peak))
      next }else{
        res <- .PICresolve(peak, pic, pval)
        pics1 <- c(pics1,list(res$pic))
        peaks1 <- c(peaks1,list(res$peak))
      }
  }

  pics$pics <- pics1
  pics$peaks <- peaks1

  return(pics)
}

.PICfit <- function(peak,pic,iter){
  # define sub-functions
  Gaussian <- function(x,position,width){
    exp(-((x-position)/(0.6005612*width))^2);
  }

  fitGaussian <- function(lambda,x,y){
    A <- matrix(0,length(x),round(length(lambda)/2))
    for(j in 1:(length(lambda)/2))
    {
      position <-lambda[2*j-1] ; width <- lambda[2*j];
      A[ ,j] <- Gaussian(x,lambda[2*j-1],lambda[2*j]);
    }
    lf=lsfit(A,y)
    PEAKHEIGHTS=abs(lf$coef[-1])
    e=y-A%*%PEAKHEIGHTS
    return(base::norm(as.matrix(e),'2'))
  }

  fitness <- function(lambda,x,y){
    -fitGaussian(lambda,x,y)
  }

  peakrange <- function(position,width,pos,wid){
    lambda <- cbind(position,width);
    NumPeaks <- nrow(lambda);
    lb <- as.vector(lambda);
    ub <- as.vector(lambda);
    for(i in 1:NumPeaks){
      ub[2*i-1] <- lambda[i,1]+pos*lambda[i,1];
      lb[2*i-1] <- lambda[i,1]-pos*lambda[i,1];
      ub[2*i] <- lambda[i,2]+wid*lambda[i,2];
      lb[2*i] <- lambda[i,2]-wid*lambda[i,2];
    }
    output <- list(ub=ub,lb=lb)
    return(output)
  }

  PEAKHEIGHTS <- function(fitresults,xa,y){
    # The positionwidth is the results of the initial estimation or GAs for position or width.
    NumPeaks <- length(fitresults)/2;
    lambda <- matrix(fitresults,nrow=NumPeaks,byrow=TRUE)
    A <- matrix(0,length(xa),round(length(lambda)/2))
    for(j in 1:(length(lambda)/2))
    {
      A[ ,j] <- Gaussian(xa,lambda[j,1],lambda[j,2]);
    }
    lf=lsfit(A,y)
    PEAKHEIGHTS=abs(lf$coef[-1])
    return(matrix(PEAKHEIGHTS))
  }

  # end of sub-functions

  hh <- 0.5*(pic[c(peak$peakIndex),2])
  ranges <- whichAsIRanges(pic[,2]>min(hh))

  if (length(peak$peakIndex)>1){
    mids <- round(diff(peak$peakIndex)/2) + peak$peakIndex[1:(length(peak$peakIndex)-1)]
    starts <- c(min(start(ranges)),mids)
    ends <- c(mids,max(end(ranges)))

    widths <- (ends-starts)*1.7
    heights <- peak$signals
    positions <- peak$peakIndex

    lambda <- as.vector(t(cbind(positions,widths)));

    lb <- peakrange(positions,widths,0.05,0.9)$lb
    ub <- peakrange(positions,widths,0.05,0.9)$ub

    GA <- ga(type = "real-valued", fitness = fitness,
             x = 1:nrow(pic), y = pic[,2], min = c(lb), max = c(ub),
             popSize=300, maxiter=iter)

    fitresults <- matrix(GA@solution, nrow=1);
    heights <- matrix(PEAKHEIGHTS(fitresults,1:nrow(pic),pic[,2]));
    positions <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s-1]})
    widths <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s]})

    fitpics <- lapply(1:length(positions),function(s){heights[s]*Gaussian(1:nrow(pic),positions[s],widths[s])})

  } else {
    starts <- min(start(ranges))
    ends <- max(end(ranges))

    widths <- (ends-starts+1)*1.7
    heights <- max(pic[,2])
    positions <- which.max(pic[,2])

    fit <- nls(pic[,2]~heights*Gaussian(1:nrow(pic),positions,widths),start=list(heights=heights,positions=positions,widths=widths),
               control=list(maxiter=500,tol=1e-3,minFactor=1/512,printEval=F,warnOnly=T))

    fitpics <- list(predict(fit,1:nrow(pic)))
    heights <- max(fitpics[[1]])
    positions <- which.max(fitpics[[1]])
    widths <- summary(fit)$parameters['widths','Estimate']
  }
  peak$peakIndex <- round(positions)
  peak$width <- round(widths)
  peak$signals <- heights
  return(list(peak=peak, fitpics=fitpics))
}

PICfit <- function(pics, iter=50){
  peaks1 <- list()
  pics1 <- list()
  for (s in 1:length(pics$pics)) {
    peak <- pics$peaks[[s]]
    pic <- pics$pics[[s]]
    res <- .PICfit(peak, pic, iter)
    if (length(res$peak$peakIndex)>1){
      for(i in 1:length(res$peak$peakIndex)){
        peak1 <- list(peakIndex=res$peak$peakIndex[i],
                      snr=res$peak$snr[i],
                      signals=res$peak$signals[i],
                      peakScale=res$peak$peakScale[i],
                      width=res$peak$width[i])
        pic1 <- cbind(pic[,1],res$fitpics[[i]])
        mz1 <- rep(NA, nrow(pic1))
        st <- max(1,round((res$peak$peakIndex[i-1]+res$peak$peakIndex[i])/2), na.rm=TRUE)
        ed <- min(nrow(pic),round((res$peak$peakIndex[i]+res$peak$peakIndex[i+1])/2), na.rm=TRUE)
        mz1[st:ed] <- pic[st:ed,3]
        pic1 <- cbind(pic1,mz1)
        colnames(pic1) <- NULL

        pics1 <- c(pics1, list(pic1))
        peaks1 <- c(peaks1, list(peak1))
      }
    } else {
      pic[,2] <- unlist(res$fitpics)
      pics1 <- c(pics1, list(pic))
      peaks1 <- c(peaks1, list(res$peak))
    }
  }
  pics$pics <- pics1
  pics$peaks <- peaks1

  return(pics)
}
