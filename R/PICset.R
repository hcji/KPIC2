readfiles <- function(files){
  filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                   "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
  filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  info <- file.info(files)
  listed <- list.files(files[info$isdir], pattern = filepattern,
                       recursive = TRUE, full.names = TRUE)
  return(listed)
}

rtequal <- function(rt0,pics){
  rt <- pics$scantime
  for (i in 1:length(pics$pics)){
    pic <- pics$pics[[i]]
    rt.i <- rt[pic[,1]]
    rt0.i <- rt0[pic[,1]]
    int.i <- approx(rt.i,pic[,2],rt0.i, rule=2)$y
    int.i[is.na(int.i)] <- 0
    pic[,2] <- int.i
    pics$pics[[i]] <- pic
  }
  pics$scantime <- rt0
  return(pics)
}

PICset <- function(files, level, mztol=0.1, gap=3, width=5, min_snr=4, equal=TRUE, export=FALSE, par=TRUE, ...){
  library(parallel)
  path <- readfiles(files)

  if (par){
    cl <- makeCluster(getOption("cl.cores", detectCores(logical = FALSE)))
    res <- parLapply(cl, path,function(filename){
      raw <- LoadData(filename)
      pics <- getPIC(raw, level, mztol, gap, width, min_snr, export)
      return(pics)
    })
    stopCluster(cl)
  } else {
    res <- lapply(path,function(filename){
      raw <- LoadData(filename)
      pics <- getPIC(raw, level, mztol, gap, width, min_snr, export)
      return(pics)
    })
  }

  if (equal){
    res <- lapply(res, function(pics){
      rtequal(res[[1]]$scantime,pics)
    })}

  return(res)
}

PICset.kmeans <- function(files, level, mztol=0.1, gap=3, width=c(5,60), min_snr=4, alpha=0.3, equal=TRUE, export=FALSE, par=TRUE, ...){
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = FALSE)))
  path <- readfiles(files)

  if (par){
    res <- parLapply(cl, path,function(filename){
      raw <- LoadData(filename)
      pics <- getPIC.kmeans(raw, level, mztol, gap, width, alpha, min_snr, export)
      return(pics)
    })
    stopCluster(cl)
  } else {
    res <- lapply(path,function(filename){
      raw <- LoadData(filename)
      pics <- getPIC.kmeans(raw, level, mztol, gap, width, alpha, min_snr, export)
      return(pics)
    })
  }

  if (equal){
    res <- lapply(res, function(pics){
      rtequal(res[[1]]$scantime,pics)
    })}

  return(res)
}

PICset.split <- function(picset, par=FALSE) {
  library(parallel)
  if (!par){
    for (i in 1:length(picset)){
      picset[[i]] <- PICsplit(picset[[i]])
    }
  }else{
    cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
    picset <- parLapply(cl,picset,function(pics){
      PICsplit(pics)
    })
    stopCluster(cl)}
  return(picset)
}

PICset.resolve <- function(picset, pval=0.01, par=FALSE) {
  library(parallel)
  if (!par){
    for (i in 1:length(picset)){
      picset[[i]] <- PICresolve(picset[[i]], pval)
    }
  }else{
    cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
    picset <- parLapply(cl,picset,function(pics){
      PICresolve(pics, pval)
    })
    stopCluster(cl)}
  return(picset)
}

PICset.fit <- function(picset, iter=50, par=FALSE) {
  library(parallel)
  if (!par){
    for (i in 1:length(picset)){
      picset[[i]] <- PICfit(picset[[i]], iter)
    }
  }else{
    cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
    picset <- parLapply(cl,picset,function(pics){
      PICfit(pics, iter)
    })
    stopCluster(cl)}
  return(picset)
}

PICset.getPeaks <- function(picset){
  for (i in 1:length(picset)){
    picset[[i]] <- getPeaks(picset[[i]])
  }
  return(picset)
}
