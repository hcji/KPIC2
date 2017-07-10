readfiles <- function(files, fullName=TRUE){
  filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                   "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
  filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  info <- file.info(files)
  listed <- list.files(files[info$isdir], pattern = filepattern,
                       recursive = TRUE, full.names = fullName)
  return(listed)
}

PICset <- function(files, level, mztol=0.1, gap=3, width=5, ...){
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
  path <- readfiles(files)
  res <- parLapply(cl, path,function(filename){
    raw <- LoadData(filename)
    pics <- getPIC(raw, level, mztol, gap, width, min_snr)
    pics$path=filename
    return(pics)
  })
  stopCluster(cl)
  return(res)
}

PICset.kmeans <- function(files, level, mztol=0.1, gap=3, width=c(5,60), alpha=0.3, ...){
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
  path <- readfiles(files)
  res <- parLapply(cl, path,function(filename){
    raw <- LoadData(filename)
    pics <- getPIC(raw, level, mztol, gap, width, alpha, min_snr)
    pics$path=filename
    return(pics)
  })
  stopCluster(cl)
  return(res)
}

PICset.detection <- function(picset, min_snr, level){
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
  res <- parLapply(cl,picset,function(pics){
    PICdetection(pics, min_snr, level)
  })
  stopCluster(cl)
  return(res)
}

PICset.split <- function(picset) {
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
  res <- parLapply(cl,picset,function(pics){
    PICsplit(pics)
  })
  stopCluster(cl)
  return(res)
}

PICset.resolve <- function(picset, pval=0.01) {
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = TRUE)))
  res <- parLapply(cl,picset,function(pics){
    PICresolve(pics, pval)
  })
  stopCluster(cl)
  return(res)
}

PICset.getPeaks <- function(picset){
  res <- lapply(picset,function(pics){
    getPeaks(pics)
  })
  return(res)
}
