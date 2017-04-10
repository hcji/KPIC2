pretreat <- function(r_fillpeaks,normalization='sum',scaling='auto scaling'){
  dat <- t(r_fillpeaks$data.mat)
  colnames(dat) <- paste('mz:',round(r_fillpeaks$features[,'mz'],2),';','rt:',round(r_fillpeaks$features[,'rt'],1))
  if (is.null(normalization)){
    dat <- dat
  }else if (normalization=='sum'){
    total <- rowSums(dat)
    dat <- dat / matrix(total,nrow(dat),ncol(dat))
  }else{
    dat <- dat / matrix(dat[,normalization],nrow(dat),ncol(dat))
  }
  
  if (scaling=='auto scaling'){
    dat <- scale(dat)
  }else if (scaling=='center'){
    dat <- scale(dat,center=T,scale=F)
  }
  return(dat)
}

peak.analyst.oplsda <- function(X, Y, ncomp, northo){
  library(ropls)
  fit <- opls(x, y, predI = ncomp, orthoI=northo)
  ValImp <- getVipVn(fit)
  ropls::plot(fit)
  return(list(result=fit,ValImp=ValImp))
}

peak.analyst.rf <- function(x, y, ntree=500, maxnodes=7){
  library(randomForest)
  y <- as.factor(y)
  fit <- randomForest(x,y,ntree=ntree,maxnodes=maxnodes,proximity=TRUE)
  ValImp <- importance(fit)
  print(fit)
  MDSplot(fit,y)
  return(list(result=fit,ValImp=ValImp))
}

peak.analyst.plsda <- function(x, y, ncomp){
  library(ropls)
  fit <- opls(x, y, predI = ncomp)
  ValImp <- getVipVn(fit)
  ropls::plot(fit)
  return(list(result=fit,ValImp=ValImp))
}

peak.analyst.koplsda <- function(X,Y,A,oax,opt='GS',nrcvouter=20,nrcvinner=10,permutation=50,kernelParams=c(1,30,29)){
  library(kopls)
  library(ggplot2)
  library(gtools)
  options(warn =-1)
  labels <- Y
  classes <- unique(Y)
  temp.Y <- data.frame(as.numeric(Y))
  Y <- matrix(0,length(Y),length(classes))
  for(i in 1:length(classes)){
    ind <- which(temp.Y==i)
    Y[ind,i] <- 1
  }
  model.opt <- koplsCVopt(X,Y,A=A,oax=oax,modelType='da',opt=opt,nrcvouter=nrcvouter,nrcvinner=nrcvinner,kernelParams=kernelParams)
  # score plot
  tmp <- data.frame(Tp.1=model.opt$koplsModel$Tp[[oax+1]][,1], To.1=model.opt$koplsModel$To[,1])
  tmp$species <- as.factor(data.frame(labels)[,])
  .theme2 <- theme(
    axis.line = element_line(colour = 'gray', size = .75), 
    panel.background = element_blank(), 
    plot.background = element_blank(),
    legend.background=element_rect(fill='white'),
    legend.key = element_blank()
  )
  points <- geom_point(data=tmp, aes_string(x=colnames(tmp)[1], y=colnames(tmp)[2], color='species'), size=3, alpha=.7)
  ell <- get.ellipse.coords(cbind(tmp[,1],tmp[,2]), group=tmp$species)# group visualization via 
  polygons <- geom_polygon(data=data.frame(ell$coords),aes(x=x,y=y, fill=group),linetype=2,alpha=0.5, show_guide = FALSE) 
  p <- ggplot(data=tmp,aes_string(x=colnames(tmp)[1], y=colnames(tmp)[2])) +
    geom_vline(xintercept = 0,linetype=2, size=.5, alpha=.5) + 
    geom_hline(yintercept = 0,linetype=2, size=.5, alpha=.5) +
    points + 
    .theme2 +
    polygons
  show(p)
  
  print("permutation testing")
  perm.y<-lapply(1:permutation,function(i)
  {
    apply(temp.Y[,1,drop=FALSE],2,gtools::permute)
  })
  if(ncol(temp.Y)==1){
    cor.with.y<-data.frame(correlation=abs(cor(cbind(temp.Y,do.call("cbind",perm.y))))[-1,1])
  } else {
    cor.with.y<-NULL
  }
  perm.model <- sapply(1:permutation,function(i){
    perm.Y <- matrix(0,length(perm.y[[i]]),length(classes))
    for(ii in 1:length(classes)){
      ind <- which(perm.y[[i]]==ii)
      perm.Y[ind,ii] <- 1
    }
    perm.model.i <- koplsCVopt(X,perm.Y,A=A,oax=oax,modelType='da',opt=opt,nrcvouter=nrcvouter,nrcvinner=nrcvinner,kernelParams=c(model.opt$KParamfinal,model.opt$KParamfinal,1))
    meanSens <- perm.model.i$da$sensSpec[[oax+1]]$totalResults$meanSens
    meanSpec <- perm.model.i$da$sensSpec[[oax+1]]$totalResults$meanSpec
    da.i <-  perm.model.i$da$sensSpec[[oax+1]]$totalResults
    ACC <- (da.i$TPtot+da.i$TNtot)/(da.i$TPtot+da.i$TNtot+da.i$FPtot+da.i$FNtot)
    return(c(meanSens,meanSpec,ACC))
  })
  permutation.test <- cbind(cor.with.y, matrix(unlist(perm.model),ncol=3))
  colnames(permutation.test) <- c('cor.with.y','meanSens','meanSpec','ACC')
  # output
  output <- list()
  da <- model.opt$da$sensSpec[[oax+1]]$totalResults
  output$ACC <- (da$TPtot+da$TNtot)/(da$TPtot+da$TNtot+da$FPtot+da$FNtot)
  output$model <- model.opt$koplsModel
  output$cv <- model.opt$cv
  output$da <- model.opt$da
  output$permutation <- permutation.test
  return(output)
}

ellipse.group<-function(obj,color,lwd=1,lty=1,border="#00000050",ellipse.level=.95,show.polygon=TRUE, alpha=.5)
{
  library(ellipse)
  library(splancs)
  
  #check color and add extra transparency
  color<-alpha.col(color,alpha)	
  #split objs and inputs for each group based on color (color, lwd, lty  are all mapped together)
  tmp.obj<-as.data.frame(obj)
  tmp.char<-as.data.frame(cbind(color,lwd,lty,border))
  fct<-as.factor(color)
  .obj<-split(tmp.obj,fct)
  .char<-split(tmp.char,fct)
  
  #calculate points for ellipse
  ellipse.var<-lapply(1:nlevels(fct),function(i)
  {
    tmp<-list()
    pts<-.obj[[i]]
    if(nrow(pts)<=2){pts<-matrix(c(NA,NA))}# avoid polygon error for 1D object
    m<-colMeans(pts)
    tmp$points<-tryCatch(ellipse(as.matrix(cov(pts)),centre=c(m[1],m[2]),level=ellipse.level),
                         error=function(e){NA})
    tmp$color<-unique(as.character(.char[[i]][,1]))[1] # choose single value
    tmp$lwd<-unique(as.numeric(as.character(.char[[i]][,2])))[1]
    tmp$lty<-unique(as.numeric(as.character(.char[[i]][,3])))[1]
    tmp$border<-unique(as.character(.char[[i]][,4]))[1]
    tmp
  })
  
  # get size to plot smallest last
  ellipse.size<-sapply(1:length(ellipse.var),function(i)
  {
    tryCatch(areapl(ellipse.var[[i]]$points),error=function(e){NA})
  })
  
  plot.order<-order(ellipse.size,decreasing=TRUE)
  
  #plot
  for(i in 1:length(ellipse.var))
  {
    if(!is.na(ellipse.var[[plot.order[i]]]$points))
    {
      if(show.polygon==TRUE)
      {
        polygon(unlist(ellipse.var[[plot.order[i]]]$points[,1]),unlist(ellipse.var[[plot.order[i]]]$points[,2]),
                col=as.character(ellipse.var[[plot.order[i]]]$color),border = ellipse.var[[plot.order[i]]]$border,
                lwd=ellipse.var[[plot.order[i]]]$lwd,lty=ellipse.var[[plot.order[i]]]$lty)		
      }else{
        lines(unlist(ellipse.var[[plot.order[i]]]$points[,1]),unlist(ellipse.var[[plot.order[i]]]$points[,2]),
              col=as.character(ellipse.var[[plot.order[i]]]$color),border = ellipse.var[[plot.order[i]]],
              lwd=ellipse.var[[plot.order[i]]]$lwd,lty=ellipse.var[[plot.order[i]]]$lty)	
      }				
    }
  }
}		

get.ellipse.coords<-function(obj,group=NULL, ellipse.level=.95){
  library(ellipse)
  library(splancs)
  
  fct<-if(is.null(group)) as.factor(rep(1,nrow(obj))) else factor(unlist(group))
  .obj<-split(as.data.frame(obj),fct)
  names<-c(colnames(obj),"")
  #calculate points for ellipse
  #for level of group
  ellipse.var<-lapply(1:nlevels(fct),function(i)
  {
    pts<-.obj[[i]]
    m<-colMeans(pts)
    tmp<-cbind(tryCatch(ellipse::ellipse(as.matrix(cov(pts)),centre=c(m[1],m[2]),level=ellipse.level),
                        error=function(e){matrix( c(NA,NA),nrow=1)}),rep(levels(fct)[i],nrow(pts)))
    colnames(tmp)<-NULL
    tmp				
  })
  
  #format for ggplot 2
  tmp<-do.call("rbind",ellipse.var)
  #remove errors
  tmp<-tmp[!is.na(tmp[,1]),]
  
  ellipse.size<-sapply(1:length(ellipse.var),function(i)
  {
    tryCatch(areapl(ellipse.var[[i]]),error=function(e){NA})
  })
  #avoiding issues with x/y as factors
  obj<-data.frame(matrix(as.numeric(as.matrix(tmp[,1:2])),ncol=2),tmp[,3])
  colnames(obj)<-c("x","y","group")
  #may need to maintain ordered factor levels
  obj$group<-factor(obj$group,levels=levels(fct),ordered=is.ordered(fct))
  return(list(coords=data.frame(obj), area=ellipse.size))	
}		