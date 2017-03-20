getDataMatrix <- function(r.group,std='maxo',min_samples=NULL){
  xset <- r.group$xset
  if (is.null(min_samples)) {
    min_samples <- round(0.7*max(r.group$xset@peakmat[,'sample']))
  }
  data.mat <- c()
  features <- c()
  group.inds <- c()
  for (i in 1:length(r.group$groups)) {
    a <- r.group$groups[[i]][[1]]
    if (length(a)<2) {next}
    features.i <- xset@peakmat[a,]
    samples <- unique(features.i[,'sample'])
    if (length(samples)<min_samples){next}
    this.vec <- rep(0,length(xset@path))
    for (j in 1:nrow(features.i)) {
      this.vec[features.i[j,'sample']] <- max(this.vec[features.i[j,'sample']],features.i[j,std])
    }
    features.i <- c(mean(features.i[,'mz']),
                    min(features.i[,'mz']),
                    max(features.i[,'mz']),
                    mean(features.i[,'rt_cor']),
                    max(features.i[,'rtmax'])-min(features.i[,'rtmin']))
    features <- rbind(features,features.i)
    data.mat <- rbind(data.mat,this.vec)
    group.inds <- c(group.inds,i)
  }
  colnames(features) <- c('mz','mzmin','mzmax','rt','width')
  colnames(data.mat) <- xset@sample
  features <- cbind(group.inds,features)
  return(list(features=features,data.mat=data.mat,r.group=r.group))
}

fillPeaks <- function(r.DataMatrix,expand_rt=1.5,expand_mz=150,min_snr=3,min_ridge=2,std='maxo'){
  xset <- r.DataMatrix$r.group$xset
  features <- r.DataMatrix$features
  data.mat <- r.DataMatrix$data.mat
  fill.inds <- which(data.mat<10^-6,arr.ind=TRUE)
  mzranges <- cbind(features[,'mzmin']-expand_mz*features[,'mz']/1000000,features[,'mzmin']+expand_mz*features[,'mz']/1000000)
  rtranges <- cbind(features[,'rt']-features[,'width']*expand_rt,features[,'rt']+features[,'width']*expand_rt)
  for (i in unique(fill.inds[,2])){
    path <- xset@path[i]
    data <- LoadData(path)
    missed <- findInterval(c(i-1,i),fill.inds[,2])
    missed <- fill.inds[(missed[1]+1):missed[2],1]
    for (j in missed) {
      this.bpc <- getBPC(data,mzranges[j,],rtranges[j,])
      r_peak_detection <- peaks_detection(this.bpc[,2],min_snr,min_ridge)
      if (length(r_peak_detection$peakIndex)==0){
        next
      }
      mainPeak <- which(abs(this.bpc[r_peak_detection$peakIndex,1]-features[j,'rt'])==min(abs(this.bpc[r_peak_detection$peakIndex,1]-features[j,'rt'])))
      r_widthEstimation <- widthEstimationCWT(this.bpc[,2],r_peak_detection)
      if (std=='signal'){
        signal <- r_peak_detection$signals[mainPeak]
        data.mat[j,i] <- signal
      } else if (std=='maxo'){
        scans_ind <- r_widthEstimation$peakIndexLower[mainPeak]:r_widthEstimation$peakIndexUpper[mainPeak]
        maxo <- max(this.bpc[scans_ind,2])
        data.mat[j,i] <- maxo
      } else if (std=='peak_area'){
        scans_ind <- r_widthEstimation$peakIndexLower[mainPeak]:r_widthEstimation$peakIndexUpper[mainPeak]
        peak_area <- integration(this.bpc[scans_ind,1],this.bpc[scans_ind,2])
        data.mat[j,i] <- peak_area
      }
    }
  }
  r.DataMatrix$data.mat <- data.mat
  return(r.DataMatrix)
}

getBPC <- function(data,mzrange,rtrange){
  scans <- findInterval(rtrange,data$times)
  scans <- (scans[1]+1):scans[2]
  bpc <- rep(0,length(scans))
  for (i in 1:length(scans)) {
    scan <- scans[i]
    peaks <- data$spectrum[[scan]]
    peak.ind <- findInterval(mzrange,peaks[,1])
    if (peak.ind[1]==peak.ind[2]){next}
    bpc[i] <- max(peaks[(peak.ind[1]+1):peak.ind[2],2])
  }
  rts <- data$times[scans]
  return(cbind(rts,bpc))
}