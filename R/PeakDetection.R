waveft = function(omega,scales)
{
  StpFrq <- omega[2]
  NbFrq  <- length(omega)
  SqrtNbFrq <- sqrt(NbFrq)
  cfsNORM <- sqrt(StpFrq)*SqrtNbFrq
  NbSc <- length(scales)
  wft <- matrix(0,NbSc,NbFrq)
  
  mul <- sqrt(scales/gamma(2+0.5))*cfsNORM
  for (jj in 1:NbSc){
    scapowered = (scales[jj]*omega)
    expnt = -(scapowered^2)/2
    wft[jj,] = mul[jj]*(scapowered^2)*exp(expnt)
  }
  return (wft)
}

cwtft = function(sig)
{
  if (length(sig)<10){stop('the length of signal is too short')}
  library(stats)
  options(warn =-1)
  val <- sig
  meanSIG <- mean(val)
  xx <- val-meanSIG
  n <- length(xx)
  dt <- 2
  
  omega <- (1:floor(n/2))
  omega <- omega*((2*pi)/(n*dt))
  omega <- c(0, omega,-(omega[seq(floor((n-1)/2),1,-1)]))
  f <- fft(xx)
  
  # set scales
  s0 <- dt
  ds <- 0.2
  NbSc <- floor(log2(n)/ds)
  scales <- c(1,s0*2^(-2:(NbSc-3)*ds))
  
  psift <- waveft(omega,scales)
  mat <- matrix(0,NbSc+1,length(f))
  for (ii in 1:(NbSc+1)){
    mat[ii,] <- f
  }
  cwtcfs <- t(mvfft(t(mat*psift),inverse = TRUE)/ncol(mat))
  cwtcfs <- cwtcfs[,1:n]
  result <- matrix(0,nrow(cwtcfs),ncol(cwtcfs))
  for (jj in 1:nrow(cwtcfs)){
    result[jj,] <- as.numeric(cwtcfs[jj,])
  }
  rownames(result) <- round(scales,2)
  return (result)
}

local_extreme <- function(data, comparator){
  mar <- matrix(0,nrow(data),2)
  mat <- cbind(mar,data,mar)
  result <- matrix(T,nrow(mat),ncol(mat))
  if (comparator=='max'){
    for (i in 3:(ncol(mat)-2)){
      Cr1 <- mat[,i]>mat[,i-2]&mat[,i]>mat[,i-1]
      Cr2 <- mat[,i]>mat[,i+2]&mat[,i]>mat[,i+1]
      result[,i] <- Cr1&Cr2
    }
  }else{if(comparator=='min'){
    for (i in 3:(ncol(mat)-2)){
      Cr1 <- mat[,i]>mat[,i-2]&mat[,i]>mat[,i-1]
      Cr2 <- mat[,i]>mat[,i+2]&mat[,i]>mat[,i+1]
      result[,i] <- Cr1&Cr2
    }
  }}
  return(result[,3:(ncol(mat)-2)])
}

ridges_detection <- function(cwt2d, sig){
  n_rows <- nrow(cwt2d)
  n_cols <- ncol(cwt2d)
  local_max <- local_extreme(cwt2d,'max')
  rows_init <- 1:5
  cols_small_peaks <- which(colSums(local_max[rows_init,])>0)
  ridges <- list()
  for (col in cols_small_peaks){
    best_rows <- rows_init[which(local_max[rows_init, col])]
    result <- ridge_detection(local_max, best_rows[1], col, n_rows, n_cols)
    rows <- result$rows
    cols <- result$cols
    staightness <- 1-sum(diff(cols))/length(cols)
    if (length(rows)>2&&staightness>0.2){
      ridges <- c(ridges,list(list(rows,cols)))
    }
  }
  ridges <- unique(ridges)
  return(ridges)
}

ridge_detection <- function(local_max, row_best, col, n_rows, n_cols){
  cols <- c()
  rows <- c()
  cols <- c(cols,col)
  rows <- c(rows,row_best)
  col_plus <- col
  col_minus <- col
  for (i in (1:n_rows)){
    row_plus <- row_best + i
    row_minus <- row_best - i
    segment_plus <- 1
    segment_minus <- 1
    if (row_minus>0 && segment_minus<col_minus && col_minus<n_cols-segment_minus-1){
      if (local_max[row_minus, col_minus + 1]){col_minus <- col_minus+1
      }else if(local_max[row_minus, col_minus - 1]){col_minus <- col_minus-1
      }else if(local_max[row_minus, col_minus]){col_minus <- col_minus
      }else{col_minus <- -1}
      if (col_minus != -1){
        rows <- c(row_minus,rows)
        cols <- c(col_minus,cols)}
    }
    if (row_plus<n_rows && segment_plus<col_plus && col_plus<n_cols-segment_plus-1){
      if (local_max[row_plus, col_plus + 1]){col_plus <- col_plus+1
      }else if(local_max[row_plus, col_plus - 1]){col_plus <- col_plus-1
      }else if(local_max[row_plus, col_plus]){col_plus <- col_plus
      }else{col_plus <- -1}
      if (col_plus != -1){
        rows <- c(rows,row_plus)
        cols <- c(cols,col_plus)}
    }
  }
  return(list(rows=rows,cols=cols))
}

peaks_position <- function(vec, ridges, cwt2d){
  n_cols <- ncol(cwt2d)
  negs <- cwt2d < 0
  local_minus <- local_extreme(cwt2d,'min')
  zero_crossing <- abs(diff(sign(cwt2d)))/2
  negs <- matrix(as.numeric(negs)+as.numeric(local_minus),nrow(local_minus),ncol(local_minus))
  negs <- negs>0
  negs[,c(0,n_cols-1)] <- T
  ridges_select <- list()
  peaks <- c()
  for (i in 1:length(ridges)){
    ridge <- ridges[[i]]
    temp <- c()
    for (j in 1:length(ridge[[1]])){
      temp <- c(temp,cwt2d[ridge[[1]][j], ridge[[2]][j]])
    }
    temp <- as.numeric(temp)
    inds <- which(temp>0)
    if (length(inds)>0){
      y <- ridge[[2]][inds]
      tb <- table(y)
      col <- round(mean(as.numeric(names(tb)[which.max(tb)])))
      rows <- ridge[[1]][(ridge[[2]] == col)]
      row <- rows[1]
      cols_start <- max(col-which(negs[row, 0:col][seq(col,1,-1)]),0)
      cols_end <- min(col+which(negs[row,col:n_cols]),n_cols)
      inds <- cols_start:cols_end
      peaks <- c(peaks,inds[vec[inds]==max(vec[inds])][1])
      ridges_select <- c(ridges_select,list(ridge))
    }else if(length(ridge[[2]])>2){
      cols_accurate <- ridge[[2]][1:(length(ridge[[2]])/2)]
      cols_start <- max(min(cols_accurate)-3, 0)
      cols_end <- min(max(cols_accurate)+4, n_cols-1)
      inds <- cols_start:cols_end
      if (length(inds)>0){
        peaks <- c(peaks,which(vec[inds]==max(vec[inds]))[1])
        ridges_select <- c(ridges_select,list(ridge))
      }
    }
  }
  return(list(peaks=peaks, ridges=ridges_select))
}

signal_noise_ratio <- function(cwt2d, ridges, peaks){
  n_cols <- ncol(cwt2d)
  row_one <- row_one_del <- cwt2d[1,]
  del <- which(abs(row_one) < 10e-5)
  if (length(del)>0){
    row_one_del <- row_one[-del]
  }
  t <- 3*median(abs(row_one_del-median(row_one_del)))/0.67
  row_one[row_one > t] <- t
  row_one[row_one < -t] <- -t
  noises <- rep(0,length(peaks))
  signals <- rep(0,length(peaks))
  scales <- rep(0,length(peaks))
  for (ind in 1:length(peaks)){
    val <- peaks[ind]
    hf_window <- length(ridges[[ind]][[2]])* 1
    win <- seq(max(c(val-hf_window, 0)),min(c(val+hf_window, n_cols)),1)
    noises[ind] <- as.numeric(quantile(abs(row_one[win]),0.95))
    for (i in 1:length(ridges[[ind]][[1]])){
      a <- ridges[[ind]][[1]][i]
      b <- ridges[[ind]][[2]][i]
      if (signals[ind] < cwt2d[a,b]){
        signals[ind] <- cwt2d[a,b]
        scales[ind] <- as.numeric(row.names(cwt2d)[a])
      }
    }
  }
  snr <- (signals+10^-5)/(noises+10^-5)
  return(list(snr=snr,signals=signals,scales=scales))
}

peaks_detection <- function(vec, min_snr, level=0, misspoint=0){
  if (length(vec)<10){
    return(NULL)
  }
  cwt2d <- cwtft(vec)
  ridges <- ridges_detection(cwt2d, vec)
  if (length(ridges)==0){
    return(NULL)
  }
  r_peaks <- peaks_position(vec, ridges, cwt2d)
  peaks <- r_peaks$peaks
  ridges <- r_peaks$ridges
  if (length(peaks)==0){return(NULL)}
  r_snr <- signal_noise_ratio(cwt2d, ridges, peaks)
  snr <- r_snr$snr
  scales <- r_snr$scales
  signals <- vec[peaks]
  peaks.limits <- c(misspoint+2,length(vec)-misspoint-2)
  
  # refine peaks
  limit1 <- which(snr>min_snr&signals>level)
  limit2 <- which(peaks>peaks.limits[1]&peaks<peaks.limits[2])
  hit <- intersect(limit2,limit1)
  
  peaks <- peaks[hit]
  snr <- snr[hit]
  scales <- scales[hit]
  signals <- signals[hit]
  
  return(list(peakIndex=peaks,snr=snr,signals=signals,peakScale=scales))
}