LoadData <- function(filename)
{
  library(mzR)
  splitname <- strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf")
  {
    msobj <- openMSfile(filename,backend="netCDF")
  }else{
    msobj <- openMSfile(filename)
  }

  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  whMS1 <- which(headerInfo$msLevel==1)
  peakInfo <- peakInfo[whMS1]

  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })

  peakNum <- unlist(lapply(peakInfo,nrow))
  scans <- unlist(lapply(1:length(peakNum),function(s){rep(s,peakNum[s])}))
  Mat <- do.call(rbind,peakInfo)

  scanTime <- round(headerInfo$retentionTime[whMS1],3)
  # close(msobj)

  return(list(path=filename, mzs=Mat$mz,ints=Mat$intensity,scans=scans,times=scanTime))
}

readPICs <- function(files){
  library(rjson)
  filepattern <- c("[j][s][o][n]")
  filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  info <- file.info(files)
  listed <- list.files(files[info$isdir], pattern = filepattern,
                       recursive = TRUE, full.names = TRUE)
  
  output <- list()
  for (i in 1:length(listed)){
    output[[i]] <- fromJSON(file=listed[i])
    gc()
  }
  return(output)
}
