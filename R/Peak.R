getNoise <- function(peaks, cwt2d, ridges){
  row_one <- row_one_del <- cwt2d[1,]
  del <- which(abs(row_one) < 10e-5)
  if (length(del)>0){
    row_one_del <- row_one[-del]
  }

  t <- 3*median(abs(row_one_del-median(row_one_del)))/0.67
  row_one[row_one > t] <- t
  row_one[row_one < -t] <- -t

  noises <- sapply(1:length(peaks),function(s){
    hf_win <- length(ridges$ridges_rows)
    win_s <- max(1, peaks[s] - hf_win)
    win_e <- min(ncol(cwt2d), peaks[s] + hf_win)
    return(as.numeric(quantile(abs(row_one[win_s:win_e]),0.9)))
  })
  return(noises)
}

peak_detection <- function(vec, min_snr, level=0){
  cwt2d <- cwtft(vec)
  sca <- cwt2d$scales
  cwt2d <- cwt2d$cwt2d
  ridges <- ridgesDetection(cwt2d, vec)
  if (length(ridges$ridges_rows)<1){return(NULL)}
  peaks <- peaksPosition(vec, ridges, cwt2d)
  signals <- getSignal(cwt2d, ridges, peaks)
  lens <- signals$ridge_lens
  lens[lens<0] <- 0
  scales <- sca[1+lens]
  lens <- signals$ridge_lens
  signals <- signals$signals
  peaks <- peaks+1
  noises <- getNoise(peaks, cwt2d, ridges)
  snr <- (signals+10^-5)/(noises+10^-5)
  refine <- snr>min_snr & lens>3 & vec[peaks]>level

  info <- cbind(peaks, scales, snr)
  info <- info[refine,]
  info <- unique(info)
  if (length(info)==0){return(NULL)
  } else if (length(info)>3){
    info <- info[order(info[,1]),]
    peakIndex=info[,1]; peakScale=info[,2]; snr=info[,3]; signals=vec[info[,1]]
  } else {
    peakIndex=info[1]; peakScale=info[2]; snr=info[3]; signals=vec[info[1]]
  }
  return(list(peakIndex=peakIndex, peakScale=peakScale, snr=snr, signals=signals))
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
}

PDM <- function(pic, min_snr, level, pval, iter){
  library(matrix)
  library(IRanges)
  vec <- pic[,2]
  rts <- pic[,1]
  mzs <- pic[,3]
  peaks <- peak_detection(vec, min_snr, level)
  
  vec.smooth <- WhittakerSmooth(vec, 2)
  helf.height <- vec.smooth[peaks$peakIndex]*0.5
  split.point <- sapply(1:(length(peaks$peakIndex)-1), function(s){
    peaks$peakIndex[s]+which.min(vec.smooth[peaks$peakIndex[s]:peaks$peakIndex[s+1]])
  })
  split.point <- unique(c(1, split.point, length(vec)))
  
  peak.ranges <- lapply(1:length(peaks$peakIndex),function(s){
    split.point[s]+which(vec.smooth[split.point[s]:split.point[s+1]]>helf.height[s])-1
  })
  
  peak.mzs <- lapply(peak.ranges, function(peak.range){
    mzs[peak.range]
  })
  
  pvals <- sapply(1:(length(peak.mzs)-1),function(s){
    t.test(peak.mzs[[s]], peak.mzs[[s+1]])$p.value
  })
  
  TP <- rep(TRUE, length(peaks$peakIndex))
  for (i in 1:length(pvals)){
    if (pvals[i] > pval){
      if (peaks$signals[i] > peaks$signals[i+1]){
        TP[i+1] <- FALSE
      } else {
        TP[i] <- FALSE
      }
    }
  }
  peaks$peakIndex <- peaks$peakIndex[TP]
  peaks$peakScale <- peaks$peakScale[TP]
  peaks$snr <- peaks$snr[TP]
  peaks$signals <- peaks$signals[TP]
  
  res <- .PICfit(peaks, pic, iter)
  plot.resolve(pic, res)
  return(res)
}

WhittakerSmooth <- function(y,lambda){
  M <- length(y)
  E <- sparseMatrix(i=1:M,j=1:M,x=1)
  D <- Matrix::diff(E)
  C <- chol(E+lambda*Matrix::t(D)%*%D)
  z <- solve(C,solve(t(C),y))
  return(as.numeric(z))
}

plot.resolve <- function(pic, res){
  library(plotly)
  rts <- pic[,1]
  raw.vec <- pic[,2]
  fit.pics <- do.call(rbind, res$fitpics)
  sum.vec <- colSums(fit.pics)
  
  p <- plot_ly(x = rts, y = raw.vec, type = 'scatter', mode = 'lines', name = 'raw') %>% 
    add_trace(x = rts, y = sum.vec, name = 'fitted') %>%
    layout(xaxis = list(title = 'Retention Time (s)'),
           yaxis = list (title = 'Intensity'))
  
  for (i in 1:nrow(fit.pics)) {
    name <- paste('peak ', i)
    p <- add_trace(p, x = rts, y = fit.pics[i,],  line = list(dash = 'dash'), name = name)
  }
  p
}
