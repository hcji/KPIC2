LoadData <- function(filename)
{
  library(mzR)
  splitname<-strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf")
  {
    mz.conn<-openMSfile(filename,backend="netCDF")
  }else{
    mz.conn<-openMSfile(filename)
  }
  
  masses<-NULL
  intensi<-NULL
  labels<-NULL
  index<-NULL
  spectrum<-list()
  b<-header(mz.conn)$retentionTime
  
  segs<-seq(0, length(b), by=200)
  if((length(b) %% 200) != 0) segs<-c(segs, length(b))
  
  for(n in 2:length(segs))
  {
    a<-mzR::peaks(mz.conn, scans=(segs[n-1]+1):segs[n])
    if (length(a)<1) {next}
    this.masses<-NULL
    this.intensi<-NULL
    this.labels<-NULL
    this.index<-NULL
    this.spectrum<-NULL
    for(i in 1:length(a))
    {
      this.a<-a[[i]]
      if (length(this.a)<1) {
        this.spectrum<-c(this.spectrum,list(this.a))
        next}
      if (!is.null(nrow(this.a)))
      {
        this.a<-this.a[this.a[,2]>1e-10,]
        if(is.null(nrow(this.a))) this.a<-matrix(this.a, nrow=1)
        this.masses<-c(this.masses, this.a[,1])
        this.intensi<-c(this.intensi, this.a[,2])
        this.labels<-c(this.labels, rep(segs[n-1]+i,nrow(this.a)))
        this.index<- c(this.index,1:length(this.a[,1]))
        this.spectrum<-c(this.spectrum,list(this.a))
      }else{
        b[segs[n-1]+i]<-NA
      }
    }
    
    masses<-c(masses, this.masses)
    intensi<-c(intensi, this.intensi)
    labels<-c(labels, this.labels)
    index<-c(index, this.index)
    spectrum<-c(spectrum,this.spectrum)
  }
  times<-b[!is.na(b)]
  close(mz.conn)
  Mat<-cbind(labels,index,masses,intensi)
  colnames(Mat)<-c('scans','index','mz','inte')
  
  return(list(Mat=Mat,spectrum=spectrum,times=times))
}

getPIC = function(filename,roi_ppm=30,level=500,itol=0.5,min_snr=3,min_ridge=3,fst=0.3,missp=5){
  library(Rcpp)
  library(stats)
  library(Ckmeans.1d.dp)
  library(forecast)
  # prepare output
  PICs <- list()
  Info <- NULL
  # load data
  data <- LoadData(filename)
  mat <- data$Mat
  rtlist <- data$times
  spectrum <- data$spectrum
  mat <- mat[order(mat[,'inte'],decreasing=T),]
  refs <- findInterval(c(-level),-mat[,'inte'])
  
  for (i in 1:refs) {
    ref.scan <- mat[i,'scans']
    ref.index <- mat[i,'index']
    ref.inte <- spectrum[[ref.scan]][ref.index,2]
    ref.mz <- spectrum[[ref.scan]][ref.index,1]
    if (length(ref.inte)<1){next}
    if (ref.inte<level){next}
    mzrange <- ref.mz*roi_ppm/1000000
    
    # set range of roi
    roi.scans <- c(1,length(rtlist))
    roi.mzs <- c(ref.mz-mzrange,ref.mz+mzrange)
    roi.mat <- NULL
    roi.inds <-NULL
    
    # locate roi
    b <- 0
    for (scan in ref.scan:roi.scans[1]){
      if (b>missp){break}
      s <- findInterval(roi.mzs,spectrum[[scan]][,1])
      if ((s[1]+1)>s[2]){
        b <- b+1
        next}
      s <- (s[1]+1):s[2]
      scan.inds <- cbind(rep(scan,length(s)),s)
      scan.mat <- spectrum[[scan]][s,]
      roi.mat <- rbind(roi.mat,scan.mat)
      roi.inds <- rbind(roi.inds,scan.inds)
      b <- 0
    }
    roi.scans[1] <- scan.inds[1]
    b <- 0
    for (scan in ref.scan:roi.scans[2]){
      if (b>missp){break}
      s <- findInterval(roi.mzs,spectrum[[scan]][,1])
      if ((s[1]+1)>s[2]){
        b <- b+1
        next}
      s <- (s[1]+1):s[2]
      scan.inds <- cbind(rep(scan,length(s)),s)
      scan.mat <- spectrum[[scan]][s,]
      roi.mat <- rbind(roi.mat,scan.mat)
      roi.inds <- rbind(roi.inds,scan.inds)
      b <- 0
    }
    roi.scans[2] <- scan.inds[1]
    # check if used
    del <- which(roi.mat[,2]<10^-6)
    if (length(del)>0){
      roi.mat <- roi.mat[-del,]
      roi.inds <- roi.inds[-del,]
    }
    if (length(roi.inds)<10){next}
    mzdiff <- abs(roi.mat[,1]-ref.mz)
    
    # kmeans cluster
    r_kmeans <- Ckmeans.1d.dp(mzdiff^2, k=c(1,5))
    mincenter <- min(r_kmeans$centers)
    tClu <- which(r_kmeans$centers==mincenter)
    select.ind <- which(r_kmeans$cluster==tClu)
    select.mat <- roi.mat[select.ind,]
    select.ind <- roi.inds[select.ind,]
    if (length(select.ind)<10){next}
    
    # refine by exponential smoothing forecasting
    PIC.scans <- roi.scans[1]:roi.scans[2]
    PIC.intensi <- rep(0,length(PIC.scans))
    start.point <- which(PIC.scans==ref.scan)
      # right side
    S1 <- S2 <- ref.inte
    miss.n <- 0
    for (j in start.point:length(PIC.scans)){
      if (miss.n>missp) {break}
      scan <- PIC.scans[j]
      
      # 2' exponential smoothing forecast
      tt <- length(S1)
      w <- fst
      at <- 2*S1[tt]-S2[tt]
      bt <- w/(1-w)*(S1[tt]-S2[tt])
      fore.intensi <- max(0,at+bt)
      
      this.scan <- which(select.ind[,1]==scan)
      if (length(this.scan)>0){
        intensi <- select.mat[this.scan,2]
        err.intensi <- (intensi-fore.intensi)/fore.intensi
        if (min(abs(err.intensi))<itol|min(err.intensi)<0) {
          miss.n <- 0
          num <- which(abs(err.intensi)==min(abs(err.intensi)))[1]
          this.scan <- this.scan[num]
          PIC.intensi[j] <- select.mat[this.scan,2]
          # remove the select ion from raw spectrum
          spectrum[[select.ind[this.scan,1]]][select.ind[this.scan,2],2] <- 0
        } else {
          PIC.intensi[j] <- fore.intensi
          miss.n <- miss.n+1
        }
      } else {
        PIC.intensi[j] <- fore.intensi
        miss.n <- miss.n+1
      }
      # 2' exponential smoothing
      S1 <- c(S1,w*PIC.intensi[j]+(1-w)*S1[tt])
      S2 <- c(S2,w*S1[tt+1]+(1-w)*S2[tt])
    }
      # left side
    S1 <- S2 <- ref.inte
    miss.n <- 0
    for (j in start.point:1){
      if (miss.n>missp) {break}
      scan <- PIC.scans[j]
      
      # 2' exponential smoothing forecast
      tt <- length(S1)
      w <- fst
      at <- 2*S1[tt]-S2[tt]
      bt <- w/(1-w)*(S1[tt]-S2[tt])
      fore.intensi <- max(0,at+bt)
      
      this.scan <- which(select.ind[,1]==scan)
      if (length(this.scan)>0){
        intensi <- select.mat[this.scan,2]
        err.intensi <- (intensi-fore.intensi)/fore.intensi
        if (min(abs(err.intensi))<itol|min(err.intensi)<0) {
          miss.n <- 0
          num <- which(abs(err.intensi)==min(abs(err.intensi)))[1]
          this.scan <- this.scan[num]
          PIC.intensi[j] <- select.mat[this.scan,2]
          # remove the select ion from raw spectrum
          spectrum[[select.ind[this.scan,1]]][select.ind[this.scan,2],2] <- 0
        } else {
          PIC.intensi[j] <- fore.intensi
          miss.n <- miss.n+1
        }
      } else {
        PIC.intensi[j] <- fore.intensi
        miss.n <- miss.n+1
      }
      # 2' exponential smoothing
      S1 <- c(S1,w*PIC.intensi[j]+(1-w)*S1[tt])
      S2 <- c(S2,w*S1[tt+1]+(1-w)*S2[tt])
    }
    
    # remove where intensity less than zero
    ll <- which(PIC.intensi>10^-6)
    if (ll[length(ll)]-ll[1]<10) {next}
    PIC.intensi <- PIC.intensi[ll[1]:ll[length(ll)]]
    PIC.scans <- PIC.scans[ll[1]:ll[length(ll)]]
    
    # peak detection
    r_peak_detection <- peaks_detection(PIC.intensi,min_snr,min_ridge,missp)
    mainPeak <- which(r_peak_detection$signal==max(r_peak_detection$signal))[1]
    if (length(r_peak_detection$peakIndex)==0){
      next
    }
    if (r_peak_detection$signals[mainPeak]<level){
      next
    }
    
    # collect infomation of PIC
    peak_rt <- rtlist[PIC.scans[r_peak_detection$peakIndex[mainPeak]]]
    peak_snr <- r_peak_detection$snr[mainPeak]
    peak_signals <- r_peak_detection$signals[mainPeak]
    peak_scale <- r_peak_detection$peakScale[mainPeak]
    r_widthEstimation <- widthEstimationCWT(PIC.intensi,r_peak_detection)
    scans_ind <- r_widthEstimation$peakIndexLower[mainPeak]:r_widthEstimation$peakIndexUpper[mainPeak]
    rtmin <- rtlist[PIC.scans[r_widthEstimation$peakIndexLower[mainPeak]]]
    rtmax <- rtlist[PIC.scans[r_widthEstimation$peakIndexUpper[mainPeak]]]
    r_peakArea <- integration(rtlist[PIC.scans[scans_ind]],PIC.intensi[scans_ind])
    mz_rsd <- sd(select.mat[,1])/mean(select.mat[,1])*1000000
    
    output <- c(mean(select.mat[,1]),
                min(select.mat[,1]),
                max(select.mat[,1]),
                peak_rt,
                rtmin,
                rtmax,
                max(PIC.intensi),
                peak_signals,
                r_peakArea,
                peak_snr,
                peak_scale,
                mz_rsd)
    # output PIC.i
    Info <- rbind(Info,output)
      # print(output)
    PIC.i <- cbind(rtlist[PIC.scans],PIC.intensi)
    colnames(PIC.i) <- c('rt','intensity')
    PICs <- c(PICs,list(PIC.i))
  }
  colnames(Info) <-  c("mz","mzmin","mzmax","rt","rtmin","rtmax","maxo","signal","peak_area","snr","scale","mz_rsd")
  index <- 1:nrow(Info)
  Info <- cbind(index,Info)
  return(list(Info=Info,PICs=PICs,rt=rtlist))
}

