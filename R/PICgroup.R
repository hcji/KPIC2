PICset.group <- function(picset, tolerance=c(0.01,10),weight=c(0.8,0.2), method='score', frac=0.5){
  library(dbscan)
  minSample = frac*length(picset)
  peakmat <- lapply(picset, function(pics){
    return(pics$peakinfo)
  })
  sample <- unlist(lapply(1:length(peakmat),function(s){
    return(rep(s,nrow(peakmat[[s]])))
  }))
  index <- unlist(lapply(1:length(peakmat),function(s){
    return(1:nrow(peakmat[[s]]))
  }))
  peakmat <- do.call(rbind,peakmat)
  peakmat <- cbind(sample,index,peakmat)
  peakmat <- peakmat[order(peakmat[,'mz']),]

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

    if (length(candidates)<minSample) {
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
  peakmat <- cbind(peakmat,group)
  peakmat <- peakmat[peakmat[,'group']!=-1,]
  return(list(peakmat=peakmat[,2:ncol(peakmat)], picset=picset))
}
