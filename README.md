# KPIC2
  KPIC2 is an effective platform for LC-MS based metabolomics using pure ion chromatograms, which is developed for metabolomics studies. KPIC2 can detect pure ions accurately, align PICs across samples, group PICs to annotate isotope and adduct PICs, fill missing peaks and pattern recognition. High-resolution mass spectrometers like TOF and Orbitrap are more suitable.

## Installation  

### Install Depends: 

    install.packages(c("devtools", "Rcpp", "RcppArmadillo", "mzR", "parallel", "shiny", "plotly", "data.table", "GA", "IRanges", "dbscan", "randomForest"))
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("mzR","ropls"))

### Install KPIC2:  

    library(devtools);  
    install_github("hcji/KPIC2")
		
## Usage:
### Feature detection:
  Extract pure ion chromatograms via optimized K-means clustering of ions in region of interest, and detect peaks of PICs.
  For a single sample:

    library(KPIC)
    filename <- 'E:/LC-MS data/example/1-1_Seg1Ev1.mzXML'
    pics <- getPIC.kmeans(filename, level=50000)
	pics <- PICsplit(pics)
	pics <- getPeaks(pics)
    viewPICs(pics)

  For a set of samples:

    files <- 'E:/LC-MS data/example'
    PICS <- PICset.kmeans(files, level=50000, export=F, par=T)
    PICS <- PICset.split(PICS)
    PICS <- PICset.getPeaks(PICS)
    viewPICs(PICS[[1]])


### Grouping and Alignment:
  Align PICs extracted in the last step across samples.

    groups_raw <- PICset.group(PICS, tolerance = c(0.1, 20))
    groups_align <- PICset.align(groups_raw, method='fftcc',move='loess')
    groups_align <- PICset.group(groups_align, tolerance = c(0.1, 20))
	groups_align <- PICset.align(groups_align, method='fftcc',move='direct')
    viewAlign(groups_raw, groups_align)
    

### Group combination:
  The group combination can be used for grouping isotopic and adduct features with the main features.

    groups_align <- groupCombine(groups_align, type='isotope')
    viewPseudospecturm(groups_align)

    
### Generate the peak matrix:
  Summarize extracted information into a data matrix.

    data <- getDataMatrix(groups_align)
    
### Miss peak filling:
  After grouping, there always be peak groups that do not in-clude peaks from every sample. The cause of missing peaks mean does lies on that the peak does not exit. It may be caused by undetection or misalignment. 

    data <- fillPeaks.EIBPC(data)
    
### Pattern recognition:
  Finding the difference between two class. random forest is used for the example.

    labels <- c(rep(1,6), rep(2,6))
    analyst.RF(labels, data)

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
