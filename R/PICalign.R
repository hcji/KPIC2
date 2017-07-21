PICset.align <- function(groups, method='fftcc', move='direct', span=1.5){
  peakmat <- groups$peakmat
  picset <- groups$picset
  rm(groups)

  ngroup <- max(peakmat[,'group'])
  id <- 1:nrow(peakmat)
  peakmat <- cbind(id, peakmat[order(peakmat[,'group']),])
  sp <- findInterval((0:ngroup)+0.5,peakmat[,'group'])
  freq <- round(mean(diff(picset[[1]]$scantime)),4)

  rt_mat <- matrix(NA, length(picset), ngroup)
  lag_mat <- matrix(NA, length(picset), ngroup)
  for (i in 1:(length(sp)-1)){
    gpi <- peakmat[(sp[i]+1):sp[i+1],]

    ref <- min(gpi[,'sample'])
    rr <- match(ref, gpi[,'sample'])
    sams <- gpi[,'sample']
    inds <- gpi[,'index']

    rsam <- sams[rr]
    rind <- inds[rr]
    rpic <- picset[[rsam]]$pics[[rind]]
    apics <- lapply(1:nrow(gpi),function(s){
      picset[[sams[s]]]$pics[[inds[s]]]
    })

    if (method == 'fftcc'){
      lags <- sapply(apics, function(apic){
        .align_fftcc(rpic, apic)
      })
    } else if (method == 'match'){
      lags <- round((gpi[rr,'rt']-gpi[,'rt'])/freq)
    } else {stop('wrong method')}

    for (j in 1:length(sams)){
      ss <- sams[j]
      rt_mat[ss,i] <- gpi[j,'rt']
      lag_mat[ss,i] <- lags[j]
    }

    if (move == 'direct') {
      for (k in 1:length(apics)){
        apics[[k]][,1] <- apics[[k]][,1] + lags[k]
        picset[[sams[k]]]$pics[[inds[k]]] <- apics[[k]]
        picset[[sams[k]]]$peakinfo[inds[k],'rt'] <- picset[[sams[k]]]$peakinfo[inds[k],'rt'] + lags[k]*freq
      }
      peakmat[gpi[,'id'],'rt'] <- peakmat[gpi[,'id'],'rt'] + lags*freq
    }
  }

  if (move == 'loess') {
    for (j in 1:nrow(lag_mat)){
      x = rt_mat[j,]
      y = lag_mat[j,]
      mod <- loess(y~x, span=span)
      indj <- peakmat[,'sample']==j
      peakinfoj <- peakmat[indj,]

      lag_p <- round(predict(mod, x=peakinfoj[,'rt']))
      lag_rt <- freq*lag_p

      peakinfoj[,'rt'] <- peakinfoj[,'rt']+lag_rt
      peakmat[indj,] <- peakinfoj
      for (k in 1:nrow(peakinfoj)){
        picset[[j]]$pics[[peakinfoj[k,'index']]][,1] <- picset[[j]]$pics[[peakinfoj[k,'index']]][,1] + lag_p[k]
        picset[[j]]$peakinfo[peakinfoj[k,'index'],'rt'] <- picset[[j]]$peakinfo[peakinfoj[k,'index'],'rt'] + lag_rt[k]
      }
    }
  }
  return(list(peakmat=peakmat, picset=picset))
}

.align_fftcc <- function(rpic,apic){
  min_scan <- min(rpic[,1], apic[,1])
  max_scan <- max(rpic[,1], apic[,1])

  ref <- rep(0, max_scan-min_scan+1)
  align <- ref

  ref[rpic[,1]-min_scan+1] <- ref[rpic[,1]-min_scan+1] + rpic[,2]
  align[apic[,1]-min_scan+1] <- align[apic[,1]-min_scan+1] + apic[,2]

  lag <- .fftcc(align, ref)
  return(lag)
}

.fftcc <- function(align, ref){
  M <- length(ref)
  X <- fft(ref)
  Y <- fft(align)
  R <- X*Conj(Y)
  R <- R/M
  rev <- fft(R,inverse=T)/length(rev)
  vals <- Re(rev)

  maxpos <- which.max(vals)
  maxi <- max(vals)

  if (maxi < 0.1){
    lag <- 0
  } else if (maxpos > length(vals)/2){
    lag = maxpos-length(vals)-1
  } else {
    lag <- maxpos-1
  }

  return(lag)
}

