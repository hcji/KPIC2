getDataMatrix <- function(r.group,std='maxo'){
  xset <- r.group$xset
  data.mat <- c()
  features <- c()
  group.inds <- c()
  for (i in 1:length(r.group$groups)) {
    a <- unlist(r.group$groups[[i]])
    features.i <- xset@peakmat[a,]
    samples <- unique(features.i[,'sample'])
    this.vec <- rep(0,length(xset@path))
    for (j in 1:nrow(features.i)) {
      this.vec[features.i[j,'sample']] <- max(this.vec[features.i[j,'sample']],features.i[j,std])
    }
    features.i <- c(mean(features.i[,'mz']),
                    min(features.i[,'mz']),
                    max(features.i[,'mz']),
                    mean(features.i[,'rt_cor']),
                    max(features.i[,'rtmax'])-min(features.i[,'rtmin']),
                    max(features.i[,'maxo']))
    features <- rbind(features,features.i)
    data.mat <- rbind(data.mat,this.vec)
    group.inds <- c(group.inds,i)
  }
  colnames(features) <- c('mz','mzmin','mzmax','rt','width','maxo')
  colnames(data.mat) <- xset@sample
  rownames(data.mat) <- NULL
  features <- cbind(group.inds,features)
  return(list(features=features,data.mat=data.mat,r.group=r.group))
}

s.fillpeaks <- function(vec,path,mzranges,rtranges,features,min_snr,min_ridge,std,peak){
  missed <- which(vec<10^-6)
  filled <- c()
  if (length(missed)<1) {return(NULL)}
  data <- LoadData(path)
  for (j in missed) {
    this.bpc <- getBPC(data,mzranges[j,],rtranges[j,])
    r_peak_detection <- peaks_detection(this.bpc[,2],min_snr,min_ridge)
    if (length(r_peak_detection$peakIndex)==0){
      filled <- c(filled,0)
      next
    }
    if (peak=='nearest'){
      mainPeak <- which(abs(this.bpc[r_peak_detection$peakIndex,1]-features[j,'rt'])==min(abs(this.bpc[r_peak_detection$peakIndex,1]-features[j,'rt'])))
    }
    if (peak=='highest'){
      mainPeak <- which(r_peak_detection$signal==max(r_peak_detection$signal))[1]
    }
    
    r_widthEstimation <- widthEstimationCWT(this.bpc[,2],r_peak_detection)
    
    if (std=='signal'){
      signal <- r_peak_detection$signals[mainPeak]
      filled <- c(filled,signal)
    } else if (std=='maxo'){
      scans_ind <- r_widthEstimation$peakIndexLower[mainPeak]:r_widthEstimation$peakIndexUpper[mainPeak]
      maxo <- max(this.bpc[scans_ind,2])
      filled <- c(filled,maxo)
    } else if (std=='peak_area'){
      scans_ind <- r_widthEstimation$peakIndexLower[mainPeak]:r_widthEstimation$peakIndexUpper[mainPeak]
      peak_area <- integration(this.bpc[scans_ind,1],this.bpc[scans_ind,2])
      filled <- c(filled,peak_area)
    }
  }
  return(list(missed=missed,filled=filled))
}

fillPeaks.peakfinder <- function(r.DataMatrix,tolerance=c(0.1,15),weight=c(0.7,0.2,0.1),std='maxo'){
  xset <- r.DataMatrix$r.group$xset
  features <- r.DataMatrix$features
  data.mat <- r.DataMatrix$data.mat
  inds <- which(data.mat<10^-3,arr.ind=TRUE)
  
  mzranges <- cbind(features[,'mz']-tolerance[1],features[,'mz']+tolerance[1])
  rtranges <- cbind(features[,'rt']-tolerance[2],features[,'rt']+tolerance[2])
  
  for (i in 1:nrow(inds)){
    mzrange.i <- mzranges[inds[i,1],]
    rtrange.i <- rtranges[inds[i,1],]
    candidates <- xset@PICset[[inds[i,2]]]$Info
    ids <- which(candidates[,'mz']>=mzrange.i[1]&candidates[,'mz']<=mzrange.i[2]&
                 candidates[,'rt']>=rtrange.i[1]&candidates[,'rt']<=rtrange.i[2])
    if (length(ids)<1){next}
    scores <- (1-abs(candidates[ids,'mz']-features[inds[i,1],'mz'])/tolerance[1])*weight[1] +
              (1-abs(candidates[ids,'rt']-features[inds[i,1],'rt'])/tolerance[2])*weight[2] +
              (1-abs(candidates[ids,'maxo']-features[inds[i,1],'maxo'])/features[inds[i,1],'maxo'])*weight[3]
    ids <- ids[which(scores==max(scores))[1]]
    data.mat[inds[1],inds[2]] <- candidates[ids,std]
  }
  return(list(features=features,data.mat=data.mat,r.group=r.group))
}

fillPeaks.EIBPC <- function(r.DataMatrix,tolerance=c(0.1,15),min_snr=3,min_ridge=2,std='maxo',peak='highest'){
  library(parallel)
  library(iterators)
  library(foreach)
  library(doParallel)
  library(KPIC)
  xset <- r.DataMatrix$r.group$xset
  features <- r.DataMatrix$features
  data.mat <- r.DataMatrix$data.mat
  
  mzranges <- cbind(features[,'mz']-tolerance[1],features[,'mz']+tolerance[1])
  rtranges <- cbind(features[,'rt']-tolerance[2],features[,'rt']+tolerance[2])
  
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = FALSE)))
  registerDoParallel(cl)
  
  result <- foreach(i=1:ncol(data.mat)) %dopar%
    s.fillpeaks(data.mat[,i],xset@path[i],mzranges,rtranges,features,min_snr,min_ridge,std,peak)
    
  stopCluster(cl)
  for (i in 1:ncol(data.mat)){
    for (j in 1:length(result[[i]]$missed)){
      this.missed <- result[[i]]$missed[j]
      this.filled <- result[[i]]$filled[j]
      r.DataMatrix$data.mat[this.missed,i] <- this.filled
    }
  }
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