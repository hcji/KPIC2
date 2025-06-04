PICset.group <- function(picset, tolerance=c(0.01,10),weight=c(0.8,0.2), method='score', frac=0.5){

  minSample = max(1, frac*length(picset))
  peakmat <- lapply(picset, function(pics){
    return(pics$peakinfo)
  })
  sample <- unlist(lapply(1:length(peakmat),function(s){
    return(rep(s,nrow(peakmat[[s]])))
  }))
  sample <- sample[sample != 0]
  index <- unlist(lapply(1:length(peakmat),function(s){
    return(seq_len(nrow(peakmat[[s]])))
  }))
  index <- index[index != 0]
  peakmat <- do.call(rbind,peakmat)
  peakmat <- cbind(sample,index,peakmat)
  peakmat <- peakmat[order(peakmat[,'mz']), , drop = FALSE]

  id <- 1:nrow(peakmat)
  group <- rep(0,nrow(peakmat))
  peakmat <- cbind(id,peakmat)

  inds <- order(-peakmat[,'maxo'])
  group_id <- 1

  for (i in inds){
    if (group[i]!=0) {next}
    mz_s <- peakmat[i,'mz']-tolerance[1]
    mz_e <- peakmat[i,'mz']+tolerance[1]
    rt_s <- peakmat[i,'rt']-tolerance[2]
    rt_e <- peakmat[i,'rt']+tolerance[2]

    mz_ref <- as.numeric(peakmat[i,'mz'])
    rt_ref <- as.numeric(peakmat[i,'rt'])

    candidates <- findCandidate(i, mz_s, mz_e, rt_s, rt_e, peakmat[,'mz'], peakmat[,'rt'], group)
    if (is.null(candidates)){next}

    if (length(candidates) <= minSample) {
      group[candidates] <- -1
      next
      } else {
        candidates <- peakmat[candidates,]
    }
    if (method=='score'){
      scores <- (1-abs(candidates[,'mz']-mz_ref)/tolerance[1])*weight[1]+
        (1-abs(candidates[,'rt']-rt_ref)/tolerance[2])*weight[2]
      hit <- unlist(lapply(1:length(picset),function(s){
        a <- which(candidates[,'sample']==s)
        b <- a[which.max(scores[a])]
        return(b)
      }))
      group[candidates[hit,'id']] <- group_id
      group_id <- group_id + 1
    } else if (method=='dbscan'){
      ref <- which(candidates[,'id']==i)
      scores <- cbind((candidates[,'mz']-mz_ref)/tolerance[1]*weight[1],
                      (candidates[,'rt']-rt_ref)/tolerance[2]*weight[2])
      r.dbscan <- hdbscan(scores,minPts=ceiling(0.5*length(picset)))
      ref.clu <- r.dbscan$cluster[ref]
      if (ref.clu==0 && sum(r.dbscan$cluster)!=0){
        group[candidates[,'id']] <- -1
        next
      } else {
        hits <- which(r.dbscan$cluster==ref.clu)
        if (length(hits)<minSample){stop()}
        candidates <- candidates[hits,]
        for (j in 1:length(picset)){
          hj <- which(candidates[,'sample']==j)
          if (length(hj)<1) {
            next
          }else if (length(hj)>1){
            hj <- hj[which.min(abs(scores[hj,1]))]}
          group[candidates[hj,'id']] <- group_id
        }
        group_id <- group_id + 1
      }
    }
  }

  peakmat <- as.data.table(cbind(peakmat,group))
  setkey(peakmat, group)

  if (all(group == -1)){
    counts <- integer(0)
    counts <- group.id <- integer(0);
    group.mz <- group.rt <- mean.ints <- numeric(0)
    peakmat <- peakmat[0]
    group.info <- cbind(group.id, group.rt, group.mz, mean.ints, counts)
  } else {
    splits <- lapply(1:max(group), function(s){
      peakmat[.(s),]
    })
    counts <- sapply(splits,nrow)
    group.id <- 1:length(splits)
    group.mz <- sapply(splits, function(s){
      round(mean(s$mz),2)
    })
    group.rt <- sapply(splits, function(s){
      round(mean(s$rt),2)
    })
    mean.ints <- sapply(splits, function(s){
      round(mean(s$maxo))
    })
    finl <- which(counts>=minSample)
    
    peakmat <- do.call(rbind, splits[finl])
    group.info <- cbind(group.id, group.rt, group.mz, mean.ints, counts)[finl,,drop=FALSE]
  }
  return(list(group.info=group.info, peakmat=peakmat[,-1], picset=picset))
}

groupCombine <- function(groups,min_corr=0.9,type='tailed',window=10){

  peakmat <- groups$peakmat
  picset <- groups$picset
  group.info <- groups$group.info
  rm(groups)

  peakmat <- as.data.table(peakmat)
  setkey(peakmat, group)

  cluster <- rep(0, nrow(group.info))
  corr <- rep(0, nrow(group.info))
  group.info <- data.table(1:nrow(group.info), group.info, cluster, corr)

  cluster.id <- 1
  for(i in 1:nrow(group.info)){
    if (group.info[i, cluster]!=0){next}
    candidate <- .groupCandidate(i, group.info, type, window)
    if (nrow(candidate)>1) {
      group.matrix <- .buildMatrix(candidate, peakmat, group.info, picset)
      group.corrs <- sapply(1:nrow(candidate), function(s){
        .mat.cor(group.matrix[[1]], group.matrix[[s]])
      })
      hits <- which(group.corrs >= min_corr)
    } else {
      group.corrs <- 1
      hits <- 1}
    hit.id <- unlist(candidate[hits,1])
    group.info[hit.id, 'cluster'] <- cluster.id
    group.info[hit.id, 'corr'] <- group.corrs[hits]
    cluster.id <- cluster.id + 1
  }

  return(list(peakmat=peakmat, picset=picset, group.info=group.info[,-1]))
}

getPseudospecturm <- function(groups, clu.id){
  peakmat <- groups$peakmat
  setkey(peakmat, group)
  group.info <- groups$group.info
  setkey(group.info, cluster)
  nsample <- length(groups$picset)
  rm(groups)

  this.groups <- group.info[.(clu.id)]
  setorder(this.groups, group.mz)
  ints.mat <- matrix(NA, nsample, nrow(this.groups))

  for (i in 1:nrow(this.groups)){
    group.id <- this.groups[i,'group.id']
    this.peaks <- peakmat[.(group.id)]
    sams <- this.peaks$sample
    ints.mat[sams,i] <- this.peaks$maxo
  }

  ints.mat <- ints.mat/ints.mat[,1]
  intensity <- colMeans(ints.mat, na.rm=TRUE) * 100
  res <- cbind(this.groups$group.mz, this.groups$group.rt, this.groups$counts, this.groups$corr, intensity)
  colnames(res) <- c('mz', 'rt', 'counts', 'corr', 'intensity')
  return(res)
}

.buildMatrix <- function(candidate, peakmat, group.info, picset){
  group.peaks <- list()
  group.matrix <- list()
  minScan <- Inf
  maxScan <- 0
  nsample <- length(picset)

  for (i in 1:nrow(candidate)){
    group.id <- candidate[i,'group.id']
    this.peaks <- peakmat[.(group.id)]
    group.peaks[[i]] <- this.peaks

    for (j in 1:nrow(this.peaks)){
      pa <- this.peaks$sample[j]
      pb <- this.peaks$index[j]
      pj <- picset[[pa]]$pics[[pb]][,1]
      minScan <- min(minScan,pj)
      maxScan <- max(maxScan,pj)
    }
  }

  nscan <- maxScan-minScan+1
  for (i in 1:nrow(candidate)){
    this.matrix <- matrix(0, nsample, nscan)
    this.peaks <- group.peaks[[i]]

    for (j in 1:nrow(this.peaks)){
      pa <- this.peaks$sample[j]
      pb <- this.peaks$index[j]
      pj <- picset[[pa]]$pics[[pb]]
      ps <- pj[,1]-minScan+1
      this.matrix[pa,ps] <- pj[,2]
    }
    group.matrix[[i]] <- this.matrix
  }
  return(group.matrix)
}

.groupCandidate <- function(i, group.info, type, window){
  ref.mz <- as.numeric(group.info[i, 'group.mz'])
  ref.rt <- as.numeric(group.info[i, 'group.rt'])

  min.rt <- ref.rt-0.5*window
  max.rt <- ref.rt+0.5*window

  if (type == 'tailed'){
    min.mz <- ref.mz
    max.mz <- ref.mz + 0.05
  } else if (type == 'isotope'){
    min.mz <- ref.mz
    max.mz <- ref.mz + 5*1.033 + 0.05
  } else if (type == 'all'){
    min.mz <- 0
    max.mz <- Inf
  }

  return(group.info[group.rt>=min.rt & group.rt<=max.rt][group.mz>=min.mz & group.mz<=max.mz][cluster==0])
}

.mat.cor <- function(A,B){
  a <- A-mean(A)
  b <- B-mean(B)
  return(sum(a*b)/sqrt(sum(a*a)*sum(b*b)))
}
