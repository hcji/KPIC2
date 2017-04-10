PICgroup <- function(xset,tolerance=c(0.01,10),weight=c(0.7,0.2,0.1),method='dbscan',frac=0.5){
  library(dbscan)
  groups <- list()
  features <- NULL
  min_samples <- round(frac*length(xset@PICset))
  peakmat <- xset@peakmat[order(xset@peakmat[,'mz']),]
  group_id <- rep(0,nrow(peakmat))
  peakmat <- cbind(1:nrow(peakmat),peakmat,group_id)
  inds <- order(-peakmat[,'maxo'])
  group.id <- 1
  for(i in inds){
    if (peakmat[i,'group_id']!=0) {next}
    group.i <- c()
    mzrange <- c(peakmat[i,'mz']-tolerance[1],peakmat[i,'mz']+tolerance[1])
    rtrange <- c(peakmat[i,'rt_cor']-tolerance[2],peakmat[i,'rt_cor']+tolerance[2])
    candidates <- findInterval(mzrange,peakmat[,'mz'])
    candidates <- (candidates[1]+1):candidates[2]
    candidates <- candidates[peakmat[candidates,'group_id']==0]
    candidates <- candidates[peakmat[candidates,'rt_cor']>=rtrange[1]&peakmat[candidates,'rt_cor']<=rtrange[2]]
    ref <- which(candidates==i)
    if (length(candidates)>=min_samples){
      candidates <- peakmat[candidates,]
    } else {peakmat[candidates,'group_id'] <- -1
    next}
    if (method=='score'){
      scores <- (1-abs(candidates[,'mz']-peakmat[i,'mz'])/tolerance[1])*weight[1] +
        (1-abs(candidates[,'rt_cor']-peakmat[i,'rt_cor'])/tolerance[2])*weight[2] +
        (1-abs(candidates[,'maxo']-peakmat[i,'maxo'])/peakmat[i,'maxo'])*weight[3]
      candidates <- cbind(candidates,scores)
      for(j in 1:length(xset@PICset)){
        sample.j <- candidates[candidates[,'sample']==j,]
        if (length(sample.j)==0){
          next
        } else if (!is.null(nrow(sample.j))){
          hit.j <- which(sample.j[,'scores']==max(sample.j[,'scores']))[1]
          group.i <- c(group.i,sample.j[hit.j,1])
          peakmat[sample.j[hit.j,1],'group_id'] <- group.id
        } else {
          group.i <- c(group.i,sample.j[1])
          peakmat[sample.j[1],'group_id'] <- group.id
        }
      }
    } else if (method=='dbscan'){
      scores <- data.frame(abs(candidates[,'mz']-peakmat[i,'mz'])/tolerance[1]*weight[1],
                           abs(candidates[,'rt_cor']-peakmat[i,'rt_cor'])/tolerance[2]*weight[2])
      r.dbscan <- hdbscan(scores,minPts=min_samples)
      ref.clu <- r.dbscan$cluster[ref]
      if (ref.clu==0&sum(r.dbscan$cluster)!=0){
        peakmat[i,'group_id']<- -1
        next}
      hits <- which(r.dbscan$cluster==ref.clu)
      
      for(j in 1:length(xset@PICset)){
        sample.j <- which(candidates[hits,'sample']==j)
        sample.j <- candidates[hits[sample.j],]
        if (length(sample.j)==0){
          next
        } else if (!is.null(nrow(sample.j))){
          hit.j <- which(sample.j[,'maxo']==max(sample.j[,'maxo']))[1]
          group.i <- c(group.i,sample.j[hit.j,1])
          peakmat[sample.j[hit.j,1],'group_id'] <- group.id
        } else {
          group.i <- c(group.i,sample.j[1])
          peakmat[sample.j[1],'group_id'] <- group.id
        }
      }
      if (length(group.i)!=length(unique(peakmat[group.i,'sample']))){stop()}
    }
    
    if (length(group.i)<min_samples){next}
    groups <- c(groups,list(group.i))
    feature.i <- c(mean(peakmat[group.i,'mz']),mean(peakmat[group.i,'rt_cor']),max(peakmat[group.i,'maxo']))
    features <- rbind(features,feature.i)
    group.id <- group.id+1
  }
  colnames(features) <- c('mz','rt','maxo')
  rownames(features) <- NULL
  xset@peakmat <- peakmat[,c(2:(ncol(peakmat)-1))]
  return(list(features=features,groups=groups,xset=xset))
}

groupCombine <- function(r.group,min_corr=0.9,type='tailed',window=10){
  groups <- r.group$groups
  features <- r.group$features
  xset <- r.group$xset
  group_id <- rep(0,nrow(features))
  features <- cbind(features,group_id)
  
  inds <- order(features[,'rt'])
  features <- features[inds,]
  groups <- groups[inds]
  combined.groups <- list()
  combined.features <- NULL
  
  inds <- order(-features[,'maxo'])
  group_id <- 1
  for(i in inds){
    if (features[i,'group_id']!=0){next}
    
    rtrange <- c(features[i,'rt']-window,features[i,'rt']+window)
    if (type=='tailed'){
      mzrange <- c(features[i,'mz'],features[i,'mz']+0.05)
    } else if (type=='isotope'){
      mzrange <- c(features[i,'mz'],features[i,'mz']+1.033*5+0.05)
    } else {
      mzrange <- c(0,Inf)
    }
    
    candidates <- findInterval(rtrange,features[,'rt'])
    candidates <- (candidates[1]+1):candidates[2]
    candidates <- candidates[features[candidates,'group_id']==0]
    candidates <- candidates[features[candidates,'mz']>=mzrange[1]&features[candidates,'mz']<=mzrange[2]]
    main_feature <- which(candidates==i)
    if (length(candidates)>1){
      PIC.matrix <- list()
      temp <- unlist(groups[candidates])
      this.rt <- c(min(xset@peakmat[temp,'rtmin']-0.5*window),max(xset@peakmat[temp,'rtmax'])+0.5*window)
      this.rt <- xset@rt[[1]][xset@rt[[1]]>=this.rt[1]&xset@rt[[1]]<=this.rt[2]]
      for (j in 1:length(candidates)){
        this.matrix <- matrix(0,length(xset@path),length(this.rt))
        this.peaks <- xset@peakmat[groups[[candidates[j]]],]
        for (k in 1:length(xset@path)){
          this.peak <- this.peaks[this.peaks[,'sample']==k,1]
          if (length(this.peak)>0){
            this.PIC <- xset@PICset[[k]]$PICs[[this.peak]]
          } else {next}
          temp.rt <- this.rt[this.rt>=this.PIC[1,1]&this.rt<=this.PIC[nrow(this.PIC),1]]
          this.PIC.y <- approx(this.PIC[,1],this.PIC[,2],temp.rt)$y
          this.matrix[k,which(this.rt%in%temp.rt)] <- this.PIC.y
        }
        PIC.matrix <- c(PIC.matrix,list(this.matrix))
      }
      coeffs <- rep(0,length(candidates))
      for (kk in 1:length(PIC.matrix)) {
        coeffs[kk] <- mat.cor(PIC.matrix[[main_feature]],PIC.matrix[[kk]])
      }
      hits <- candidates[coeffs>=min_corr]
      features[hits,'group_id'] <- group_id
      if (type=='tailed'){
        combined.groups <- c(combined.groups, list(groups[i]))
      } else {
        combined.groups <- c(combined.groups, list(groups[c(i,setdiff(hits,i))]))
      }
      combined.features <- rbind(combined.features,features[i,])
      group_id <- group_id + 1
    }else{
      features[i,'group_id'] <- group_id
      combined.groups <- c(combined.groups, groups[i])
      combined.features <- rbind(combined.features,features[i,])
      group_id <- group_id + 1
    }
  }
  groups <- combined.groups
  features <- combined.features
  return(list(features=features,groups=groups,xset=xset))
}

getPseudospecturm <- function(r.group,index,std='maxo') {
  peakmat <- r.group$xset@peakmat
  if (length(r.group$groups[[index]])<2) {return(NULL)}
  intensi.mat <- matrix(0,length(xset@sample),length(r.group$groups[[index]]))
  mz.mat <- matrix(NA,length(xset@sample),length(r.group$groups[[index]]))
  for (i in 1:ncol(intensi.mat)) {
    peaks.i <- r.group$groups[[index]][[i]]
    for (peak in peaks.i) {
      intensi.mat[peakmat[peak,'sample'],i] <- max(intensi.mat[peakmat[peak,'sample'],i],peakmat[peak,std])
      mz.mat[peakmat[peak,'sample'],i] <- peakmat[peak,'mz']
    }
  }
  Pseudospecturm <- c()
  for (j in 1:nrow(intensi.mat)) {
    if (prod(intensi.mat[j,])==0) {next}
    Pseudospecturm.j <- intensi.mat[j,-1]/intensi.mat[j,1]
    Pseudospecturm <- rbind(Pseudospecturm,Pseudospecturm.j)
  }
  Pseudospecturm <- cbind(rep(1,nrow(Pseudospecturm)),Pseudospecturm)
  Pseudospecturm <- colMeans(Pseudospecturm*100)
  mz.mat <- colMeans(mz.mat,na.rm = TRUE)
  Pseudospecturm <- data.frame(mz.mat,Pseudospecturm)
  Pseudospecturm <- Pseudospecturm[order(Pseudospecturm[,1]),]
  colnames(Pseudospecturm) <- c('mz','intensity')
  stem(Pseudospecturm[,1],Pseudospecturm[,2])
  return(Pseudospecturm)
}
mat.cor <- function(A,B){
  a <- A-mean(A)
  b <- B-mean(B)
  return(sum(a*b)/sqrt(sum(a*a)*sum(b*b)))
}