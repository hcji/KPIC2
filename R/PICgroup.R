PICgroup <- function(xset, beta=0.8, group.width=15, ppm=80, min_CoP=0.9){
  peakmat <- xset@peakmat
  min_samples <- round(0.7*max(peakmat[,'sample']))
  group_id <- rep(0,nrow(peakmat))
  peakmat <- cbind(peakmat,group_id)
  peakmat <- peakmat[order(peakmat[,'rt_cor']),]
  intensi_order <- order(peakmat[,'maxo'],decreasing=T)
  
  pseudospecturm <- list()
  groups <- list()
  data.mat <- c()
  group_rts <- c()
  group_ind <- 1
  for (i in 1:length(intensi_order)){
    base.feature <- intensi_order[i]
    if (peakmat[base.feature,'group_id']!=0){next}
    if (is.null(group.width)){
      FWHM <- max(0.1,0.58875*(peakmat[base.feature,'rtmax']-peakmat[base.feature,'rtmin']))
      group.width <- beta*FWHM
    }
    candi.feature <- findInterval(c(peakmat[base.feature,'rt_cor']-group.width,peakmat[base.feature,'rt_cor']+group.width),peakmat[,'rt_cor'])
    candi.feature <- (candi.feature[1]+1):candi.feature[2]
    candi.feature <- candi.feature[peakmat[candi.feature,'group_id']==0]
    candi.mz <- peakmat[candi.feature,'mz']
    candi.mz.order <- order(candi.mz)
    candi.mzdiff <- 1000000*diff(candi.mz[candi.mz.order])/candi.mz[candi.mz.order][1:(length(candi.mz.order)-1)]
    b.points <- which(candi.mzdiff>ppm)
    b.points <- candi.mz[candi.mz.order][b.points+1]
    candi.clu <- 1+findInterval(candi.mz,b.points)
    main.feature <- 1+findInterval(peakmat[base.feature,'mz'],b.points)
    feature.num <- table(candi.clu)
    hits <- which(feature.num>min_samples)
    if ((!main.feature%in%hits)|length(hits)<2) {
      peakmat[candi.feature[candi.clu%in%main.feature],'group_id'] <- group_ind
      group_ind <- group_ind+1
      groups <- c(groups,list(candi.feature[candi.clu==main.feature]))
      group_rts <- c(group_rts, peakmat[base.feature,'rt_cor'])
      next
    }
    candi.feature <- candi.feature[candi.clu%in%hits]
    candi.clu <- candi.clu[candi.clu%in%hits]
    rt.range <- c(min(peakmat[candi.feature,'rtmin']),max(peakmat[candi.feature,'rtmax']))
    rtlist <- seq(rt.range[1],rt.range[2],xset@rt[[1]][2]-xset@rt[[1]][1])
    PIC.matrix <- list()
    for (clu in hits) {
      this.PIC.matrix <- matrix(0,length(xset@path),length(rtlist))
      this.clu <- candi.feature[candi.clu%in%clu]
      this.sample <- unique(peakmat[this.clu,'sample'])
      for (sam in this.sample) {
        this.peak <- this.clu[peakmat[this.clu,'sample']==sam]
        this.peak <- this.peak[peakmat[this.peak,'maxo']==max(peakmat[this.peak,'maxo'])][1]
        this.index <- peakmat[this.peak,'index']
        this.PIC <- xset@PICset[[sam]]$PICs[[this.index]]
        this.PIC.matrix[sam,] <- approx(this.PIC[,1],this.PIC[,2],rtlist)$y
      }
      this.PIC.matrix[is.na(this.PIC.matrix)] <- 0
      PIC.matrix <- c(PIC.matrix,list(this.PIC.matrix))
    }
    
    coeffs <- rep(0,length(hits))
    reff.clu <- which(hits==main.feature)
    for (j in 1:length(PIC.matrix)) {
      coeffs[j] <- mat.cor(PIC.matrix[[reff.clu]],PIC.matrix[[j]])
    }
    select.clu <- hits[coeffs>min_CoP]
    peakmat[candi.feature[candi.clu%in%select.clu],'group_id'] <- group_ind
    group_ind <- group_ind+1
    this.groups <- list(candi.feature[candi.clu%in%main.feature])
    if (length(select.clu)>1){
      for (s.clu in setdiff(select.clu,main.feature)) {
        this.groups <- c(this.groups,list(candi.feature[candi.clu%in%s.clu]))
      }
    }
    groups <- c(groups,list(this.groups))
    group_rts <- c(group_rts, peakmat[base.feature,'rt_cor'])
  }
  xset@peakmat <- peakmat
  return(list(groups=groups,group.rts=group_rts,xset=xset))
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