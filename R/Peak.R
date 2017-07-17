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

peak_detection <- function(vec, min_snr, scales, level=0){
  cwt2d <- cwtft(vec)
  sca <- cwt2d$scales
  cwt2d <- cwt2d$cwt2d
  ridges <- ridgesDetection(cwt2d, vec)
  peaks <- peaksPosition(vec, ridges, cwt2d)
  signals <- getSignal(cwt2d, ridges, peaks)
  scales <- sca[1+signals$ridge_lens]
  signals <- signals$signals
  peaks <- peaks+1
  noises <- getNoise(peaks, cwt2d, ridges)
  snr <- (signals+10^-5)/(noises+10^-5)
  refine <- snr>min_snr & scales>4 & vec[peaks]>level

  peaks <- peaks[refine]
  scales <- scales[refine]
  snr <- snr[refine]
  info <- cbind(peaks, scales, snr)
  info <- unique(info)
  return(list(peakIndex=info[,1], scales=info[,2], snr=info[,3]))
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
}
