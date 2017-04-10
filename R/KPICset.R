setClass("KPICSet", representation(peakmat="matrix",
                                   rt="list",
                                   n_features="numeric",
                                   sample="vector",
                                   PICset="list",
                                   path="character"))

KPICset <-
  function(files,roi_range=0.1,level=500,itol=c(-0.5,0.3),min_snr=3,peakwidth=c(5,60),min_ridge=3,fst=0.3,missp=5,cluster.ref="square",eval=TRUE){
    library(parallel)
    library(iterators)
    library(foreach)
    library(doParallel)
    library(Ckmeans.1d.dp)
    library(compiler)
    library(KPIC)
    output <- new("KPICSet")
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    samples <- list.files(files[info$isdir], pattern = filepattern,
                          recursive = FALSE, full.names = FALSE)
    
    cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
    registerDoParallel(cl)
    
    result <- foreach(i=1:length(listed)) %dopar%
      getPIC(listed[i],roi_range,level,itol,min_snr,peakwidth,min_ridge,fst,missp,cluster.ref)
    
    if (eval){
      peakmat <- foreach(i=1:length(listed)) %dopar%
        PICeval(result[[i]])
      peakmat <- do.call(rbind, peakmat)
    }else{
      peakmat <- c()
      for (i in 1:length(listed)){
        peakmat <- rbind(peakmat,result[[i]]$Info)
      }
    }
    sample <- foreach(i=1:length(listed)) %do%
      rep(i,nrow(result[[i]]$Info))
    n_features <- foreach(i=1:length(listed)) %do%
      nrow(result[[i]]$Info)
    n_features <- unlist(n_features)
    sample <- unlist(sample)
    peakmat <- cbind(peakmat,sample)
    
    output@peakmat <- peakmat
    output@n_features <- n_features
    output@sample <- samples
    output@path <- listed
    output@PICset <- result
    output@rt <- foreach(i=1:length(listed)) %do%
      result[[i]]$rt
    stopCluster(cl)
    return(output)
  }

PIAlign <- function(xset,mzSegSize=.5,shift=10,lambda=1.5,ref=NA){
  library(KPIC)
  library(Matrix)
  lambda <- max(0,lambda)
  path <- xset@path
  if (is.na(ref)){
    TICs <- getTIC(path,Ref=TRUE)
    ref <- TICs$ref
  }
  peakmat <- xset@peakmat
  shift <- round(shift/(xset@rt[[1]][2]-xset@rt[[1]][1]))
  mzlist1 <- seq((min(xset@peakmat[,'mz'])-0.5*mzSegSize),
                 (max(xset@peakmat[,'mz'])+0.5*mzSegSize),mzSegSize)
  mzlist2 <- mzlist1+0.5*mzSegSize
  peakmat[,'rt'] <- rt_cor <- round(peakmat[,'rt'],2)
  peakmat <- cbind(peakmat,rt_cor)
  xset@peakmat <- peakmat
  for (i in 1:(length(mzlist1)-1)){
    xset <- seg_rtcor(xset,ref,mzlist1[i],mzlist1[i+1],shift=shift,lambda=lambda)
    xset <- seg_rtcor(xset,ref,mzlist2[i],mzlist2[i+1],shift=shift,lambda=lambda)
    if (i%%20==0){cat(round(i/(length(mzlist1)-1)*100),'%',' is done','\n')}
  }
  return(xset)
}

getTIC <- function(path,Ref=TRUE)
{
  library(KPIC)
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical=FALSE)))
  TICs <- parLapply(cl,path,function(path){
    data <- KPIC::LoadData(path)
    inte <- rep(0,length(data$times))
    for (j in 1:length(data$times)){
      inte[j] <- sum(data$spectrum[[j]][,2])
    }
    return(inte)
  })
  stopCluster(cl)
  ref <- NA
  if (Ref){
    scans <- 10^6
    corrs <- matrix(0,length(TICs),length(TICs))
    corrs1 <- rep(0,length(TICs))
    for (i in 1:length(TICs)){
      for (j in 1:length(TICs)){
        temp <- min(length(TICs[[i]]),length(TICs[[j]]))
        corrs[i,j] <- cor(TICs[[i]][1:temp],TICs[[j]][1:temp])
        scans <- min(scans,temp)
      }
      corrs1[i] <- mean(corrs[i,])
    }
    ref <- which(corrs1==max(corrs1))
  }
  return(list(TICs=TICs,ref=ref))
}