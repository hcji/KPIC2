getDataMatrix <- function(groups, std='maxo'){
  library(data.table)
  group.info <- as.data.table(groups$group.info)
  if (!'cluster'%in%colnames(group.info)){
    group.info[,cluster:=1:nrow(group.info)]
  }
  setkey(group.info, cluster)
  peakmat <- as.data.table(groups$peakmat)
  setkey(peakmat, group)
  nsam <- length(groups$picset)
  nvar <- max(group.info$cluster)

  data.mat <- matrix(0, nsam+6, nvar)
  sample.name <- NULL
  for (i in 1:nsam){
    path <- groups$picset[[i]]$path
    str1 <- unlist(strsplit(path, '/'))
    str1 <- str1[length(str1)]
    str2 <- strsplit(str1,'\\.')
    str2 <- str2[[1]][1]
    sample.name <- c(sample.name, str2)
  }

  rownames(data.mat) <- c('rt','mz','rtmin','rtmax','mzmin','mzmax',sample.name)
  info <- NULL
  for (i in 1:nvar){
    this.var <- group.info[.(i)]
    this.group <- this.var[which.max(this.var$mean.ints)]$group.id
    peaks <- as.data.frame(peakmat[.(this.group)])

    rt <- round(mean(peaks$rt),2)
    rtmin <- min(peaks$rtmin)
    rtmax <- max(peaks$rtmax)
    mz <- round(mean(peaks$mz),4)
    mzmin <- min(peaks$mzmin)
    mzmax <- max(peaks$mzmax)

    a <- peaks[,'sample']
    b <- round(peaks[,std])
    data.mat[1:6,i] <- c(rt,mz,rtmin,rtmax,mzmin,mzmax)
    data.mat[a+6,i] <- b
  }

  groups$data.mat <- data.mat
  return(groups)
}

fillPeaks.EIBPC <- function(groups, extand_mz=20, extand_rt=5, min_snr=3, std='maxo'){
  nsam <- length(groups$picset)
  nvar <- ncol(groups$data.mat)
  data.mat <- groups$data.mat

  for (i in 1:nsam){
    blank <- which(data.mat[i+6,]==0)
    if (length(blank)<1){next}
    path <- groups$picset[[i]]$path
    cat(paste('filling peaks of', path))
    raw <- LoadData(path)
    for (j in blank){
      rtrange <- c(data.mat['rtmin',j]-extand_rt, data.mat['rtmax',j]+extand_rt)
      mzrange <- c(data.mat['mzmin',j]*(1-extand_mz/10^6), data.mat['mzmax',j]*(1+extand_mz/10^6))
      eibpc <- .EIBPC(raw, mzrange, rtrange)
      ifpeak <- peak_detection(eibpc[,2], min_snr)
      if (length(ifpeak$peakIndex)>0){
        if (std=='maxo'){this.int <- max(eibpc[,2])
        }else if (std=='area') {this.int <- integration(raw$times[eibpc[,1]],eibpc[,2])}
        data.mat[i+6,j] <- this.int
      }
    }
  }
  groups$data.mat <- data.mat
  return(groups)
}

.EIBPC <- function(raw, mzrange, rtrange){
  library(data.table)
  scans <- findInterval(rtrange, raw$times)
  inds <- findInterval(scans, raw$scans)
  inds <- inds[1]:inds[2]
  rtroi <- data.table(round(10^4*raw$mzs[inds]), raw$ints[inds], raw$scans[inds])
  colnames(rtroi) <- c('mzs','ints','scans')
  setkey(rtroi, mzs)

  mzrange <- round(10^4*mzrange)
  roi <- rtroi[.(mzrange[1]:mzrange[2]),nomatch=FALSE]
  setkey(roi,scans)

  scan.bpc <- roi$scans[1]:roi$scans[length(roi$scans)]
  ints.bpc <- rep(0, length(scan.bpc))
  for (i in 1:length(scan.bpc)){
    this.scan <- scan.bpc[i]
    this.ints <- max(0, roi[.(this.scan)]$ints, na.rm=TRUE)
    ints.bpc[i] <- this.ints
  }

  res <- cbind(scan.bpc, ints.bpc)
  return(res)
}
