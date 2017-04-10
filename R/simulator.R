simulator <-  function(ground_truth,noise_file,output_file,width=c(10,50),height=c(10,10000),snratio=Inf,noise.level=0,ppm=0){
  library(readr)
  library(caMassClass)
  raw <- as.data.frame(read_csv(ground_truth))
  widths <- runif(nrow(raw),width[1],width[2])
  heights <- runif(nrow(raw),height[1],height[2])
  rawmzxml <- LoadData(noise_file)
  scantime <- rawmzxml$times
  mat <- c()
  for (i in 1:nrow(raw)) {
    scan.i <- which(abs(scantime-raw[i,'retentiontime'])==min(abs(scantime-raw[i,'retentiontime'])))
    sigma.i <- widths[i]
    scans <- max(1,round(scan.i-widths[i]/0.25)):min(length(scantime),round(scan.i+widths[i]/0.25))
    intensities <- heights[i]*exp(-(scans-scan.i)^2/(sigma.i^2))
    intensities <- intensities + 1/snratio*rnorm(length(intensities))*intensities
    mzs <- rep(raw[i,'mzs'],length(scans)) 
    mzs <- mzs + runif(length(scans),-1,1)*raw[i,'mzs']*ppm/1000000
    
    mat.i <- cbind(scans,mzs,intensities)
    mat <- rbind(mat,mat.i)
  }
  
  noise.mat <- rawmzxml$Mat[,c(1,3,4)]
  noise.mat[,3] <- abs(rnorm(nrow(noise.mat)))*noise.level
  
  fin.mat <- rbind(mat,noise.mat)
  fin.mat <- as.data.frame(fin.mat)
  fin.mat <- fin.mat[order(fin.mat$scans),]
  
  mzXML = read.mzXML(noise_file)
  for (j in 1:length(mzXML$scan)) {
    inds <- findInterval(c(j-1,j),fin.mat$scans)
    inds <- (inds[1]+1):inds[2]
    temps <- fin.mat[inds,c(2,3)]
    temps <- temps[order(temps$mzs),]
    mass <- temps$mzs
    peaks <- temps$intensities
    mzXML$scan[[j]]$mass <- mass
    mzXML$scan[[j]]$peaks <- peaks
  }
  
  write.mzXML(mzXML, output_file, precision='32') 
  cat('finish')
}