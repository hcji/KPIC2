LoadData <- function(filename)
{
  library(mzR)
  splitname <- strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf")
  {
    msobj <- openMSfile(filename,backend="netCDF")
  }else{
    msobj <- openMSfile(filename)
  }
  
  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  whMS1 <- which(headerInfo$msLevel==1)
  peakInfo <- peakInfo[whMS1]
  
  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })
  
  peakNum <- unlist(lapply(peakInfo,nrow))
  index <- unlist(lapply(peakNum,function(PN){1:PN}))
  scans <- unlist(lapply(1:length(peakNum),function(s){rep(s,peakNum[s])}))
  Mat <- cbind(scans,index, do.call(rbind,peakInfo))
  
  scanTime <- headerInfo$retentionTime[whMS1]
  close(msobj)

  return(list(Mat=Mat,spectrum=peakInfo,times=scanTime))
}

locateROI <- function(spectrum,ref.scan,roi.scans,roi.mzs,missp){
  roi.mat <- list()
  roi.inds <- list()
  roi.scans <- roi.scans[1]:roi.scans[2]
  p <- match(ref.scan,roi.scans)
  # left direction
  b <- 0
  for (s in p:1){
    if (b>missp){break}
    scan <- roi.scans[s]
    inc <- findInterval(roi.mzs,spectrum[[scan]][,1])
    if (inc[2]<=inc[1]){
      b <- b+1
      next}
    inc <- (inc[1]+1):inc[2]
    inc <- inc[spectrum[[scan]][inc,2]>0]
    if (length(inc)<1) {
      b <- b+1
      next
    } # check if used
    scan.inds <- cbind(rep(scan,length(inc)),inc)
    scan.mat <- spectrum[[scan]][inc,]
    roi.mat[[s]] <- scan.mat
    roi.inds[[s]] <- scan.inds
    b <- 0
  }

  # right direction
  b <- 0
  for (s in p:length(roi.scans)){
    if (b>missp){break}
    scan <- roi.scans[s]
    inc <- findInterval(roi.mzs,spectrum[[scan]][,1])
    if (inc[2]<=inc[1]){
      b <- b+1
      next}
    inc <- (inc[1]+1):inc[2]
    inc <- inc[spectrum[[scan]][inc,2]>0]
    if (length(inc)<1) {
      b <- b+1
      next
    } # check if used
    scan.inds <- cbind(rep(scan,length(inc)),inc)
    scan.mat <- spectrum[[scan]][inc,]
    roi.mat[[s]] <- scan.mat
    roi.inds[[s]] <- scan.inds
    b <- 0
  }
  
  roi.mat <- do.call(rbind,roi.mat)
  roi.inds <- do.call(rbind,roi.inds)
  roi.scans <- c(roi.inds[1,1],roi.inds[nrow(roi.inds),1])
  return(list(roi.scans=roi.scans, roi.inds=roi.inds, roi.mat=roi.mat))
}

ionRefine <- function(PIC.scans,ref.scan,ref.inte,select.ind,select.mat,spectrum,fst){
  PIC.intensi <- rep(0,length(PIC.scans))
  start.point <- which(PIC.scans==ref.scan)
  # right side
  S1 <- S2 <- ref.inte
  for (j in start.point:length(PIC.scans)){
    scan <- PIC.scans[j]
    # 2' exponential smoothing forecast
    tt <- length(S1)
    w <- fst
    at <- 2*S1[tt]-S2[tt]
    bt <- w/(1-w)*(S1[tt]-S2[tt])
    fore.intensi <- max(0,at+bt)
    
    this.scan <- which(select.ind[,1]==scan)
    if (length(this.scan)>1){
      intensi <- select.mat[this.scan,2]
      p <- this.scan[which.min(abs(intensi-fore.intensi))]
      PIC.intensi[j] <- select.mat[p,2]
      spectrum[[select.ind[p,1]]][select.ind[p,2],2] <- 0
    } else if(length(this.scan)==1){
      p <- this.scan
      PIC.intensi[j] <- select.mat[p,2]
      spectrum[[select.ind[p,1]]][select.ind[p,2],2] <- 0
    }else{
      PIC.intensi[j] <- fore.intensi
    }
    # 2' exponential smoothing
    S1 <- c(S1,w*PIC.intensi[j]+(1-w)*S1[tt])
    S2 <- c(S2,w*S1[tt+1]+(1-w)*S2[tt])
  }
  
  # right side
  S1 <- S2 <- ref.inte
  for (j in start.point:1){
    scan <- PIC.scans[j]
    # 2' exponential smoothing forecast
    tt <- length(S1)
    w <- fst
    at <- 2*S1[tt]-S2[tt]
    bt <- w/(1-w)*(S1[tt]-S2[tt])
    fore.intensi <- max(0,at+bt)
    
    this.scan <- which(select.ind[,1]==scan)
    if (length(this.scan)>1){
      intensi <- select.mat[this.scan,2]
      p <- this.scan[which.min(abs(intensi-fore.intensi))]
      PIC.intensi[j] <- select.mat[p,2]
      spectrum[[select.ind[p,1]]][select.ind[p,2],2] <- 0
    } else if(length(this.scan)==1){
      p <- this.scan
      PIC.intensi[j] <- select.mat[p,2]
      spectrum[[select.ind[p,1]]][select.ind[p,2],2] <- 0
    }else{
      PIC.intensi[j] <- fore.intensi
    }
    # 2' exponential smoothing
    S1 <- c(S1,w*PIC.intensi[j]+(1-w)*S1[tt])
    S2 <- c(S2,w*S1[tt+1]+(1-w)*S2[tt])
  }
  
  return(list(spectrum=spectrum,PIC.scans=PIC.scans,PIC.intensi=PIC.intensi))
}

getPIC = function(filename,roi_range=0.1,level=500,min_snr=3,peakwidth=c(5,60),fst=0.3,missp=5){
  library(Ckmeans.1d.dp)
  library(data.table)
  # prepare output
  PICs <- list()
  Info <- NULL
  # load data
  data <- LoadData(filename)
  mat <- data$Mat
  rtlist <- data$times
  spectrum <- data$spectrum
  rm(data)
  
  min_width <- round(peakwidth[1]/(rtlist[2]-rtlist[1]))
  max_width <- round(peakwidth[2]/(rtlist[2]-rtlist[1])/2)
  mzrange <- roi_range/2
  
  # set seeds
  mat <- mat[mat[,'intensity']>=level,]
  mat <- mat[order(mat[,'intensity']),]

  for (i in 1:nrow(mat)) {
    ref.scan <- as.numeric(mat[i,'scans'])
    ref.index <- as.numeric(mat[i,'index'])
    ref.inte <- spectrum[[ref.scan]][ref.index,2]
    ref.mz <- spectrum[[ref.scan]][ref.index,1]
    if (length(ref.inte)<1){next}
    if (ref.inte<level){next}

    # set range of roi
    roi.scans <- c(max(1,ref.scan-max_width),min(length(rtlist),ref.scan+max_width))
    roi.mzs <- c(ref.mz-mzrange,ref.mz+mzrange)
    roi.mat <- NULL
    roi.inds <-NULL
    
    # locate roi
    roi <- locateROI(spectrum,ref.scan,roi.scans,roi.mzs,missp)
    roi.scans <- roi$roi.scans
    roi.inds <- roi$roi.inds
    roi.mat <- roi$roi.mat
    rm(roi)
    
    # check roi length
    if (roi.scans[2]-roi.scans[1]<min_width){
      spectrum[[ref.scan]][ref.index,2] <- 0
      next}
    
    # calculate m/z difference
    mzdiff <- (roi.mat[,1]-ref.mz)^2
    
    # kmeans cluster
    r_kmeans <- Ckmeans.1d.dp(mzdiff, k=c(1,5))
    mincenter <- min(r_kmeans$centers)
    tClu <- which(r_kmeans$centers==mincenter)
    sel <- which(r_kmeans$cluster==tClu)
    if (length(sel)<min_width){next}
    select.mat <- roi.mat[sel,]
    select.ind <- roi.inds[sel,]
    
    # refine by exponential smoothing forecasting
    PIC.scans <- roi.scans[1]:roi.scans[2]
    pic <- ionRefine(PIC.scans,ref.scan,ref.inte,select.ind,select.mat,spectrum,fst)
    spectrum <- pic$spectrum
    PIC.scans <- pic$PIC.scans
    PIC.intensi <- pic$PIC.intensi
    rm(pic)
    
    # peak detection
    r_peak_detection <- peaks_detection(PIC.intensi,min_snr,level,missp)
    mainPeak <- which.max(r_peak_detection$signal)
    if (length(r_peak_detection$peakIndex)==0){
      next
    }
    if (r_peak_detection$signals[mainPeak]<level){
      next
    }
    
    # collect infomation of PIC
    peak_rt <- rtlist[PIC.scans[r_peak_detection$peakIndex[mainPeak]]]
    peak_snr <- r_peak_detection$snr[mainPeak]
    peak_scale <- r_peak_detection$peakScale[mainPeak]
    rtmin <- rtlist[PIC.scans[1]]
    rtmax <- rtlist[PIC.scans[length(PIC.scans)]]
    mz_rsd <- sd(select.mat[,1])/mean(select.mat[,1])*1000000
    
    output <- c(ref.mz,
                min(select.mat[,1]),
                max(select.mat[,1]),
                peak_rt,
                rtmin,
                rtmax,
                max(PIC.intensi),
                peak_snr,
                peak_scale,
                mz_rsd)
    # output PIC.i
    Info <- rbind(Info,output)
    PIC.i <- cbind(rtlist[PIC.scans],PIC.intensi)
    colnames(PIC.i) <- c('rt','intensity')
    PICs <- c(PICs,list(PIC.i))
  }
  colnames(Info) <-  c("mz","mzmin","mzmax","rt","rtmin","rtmax","maxo","snr","scale","rsd")
  index <- 1:nrow(Info)
  Info <- cbind(index,Info)
  gc()
  return(list(Info=Info,PICs=PICs,rt=rtlist))
}

