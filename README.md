# KPIC2
  KPIC2 is an effective platform for LC-MS based metabolomics using pure ion chromatograms, which is developed for metabolomics studies. KPIC2 can detect pure ions accurately, align PICs across samples, group PICs to annotate isotope and adduct PICs, fill missing peaks and pattern recognition. High-resolution mass spectrometers like TOF and Orbitrap are more suitable.

## Installation  

### Install Depends: 

    install.packages(c("BiocManager", "devtools", "Ckmeans.1d.dp", "Rcpp", "RcppArmadillo", "mzR", "parallel", "shiny", "plotly", "data.table", "GA", "IRanges", "dbscan", "randomForest"))
    BiocManager::install(c("mzR","ropls"))

### Install KPIC2:  

#### Released version (Suggested)
Download the source package at [url](https://github.com/hcji/KPIC2/releases) and install the package locally.

#### Development version
Note: Development may include bugs/error and functions not ready.

    library(devtools);  
    install_github("hcji/KPIC2")
	
## Usage:
### Feature detection:
  Extract pure ion chromatograms via optimized K-means clustering of ions in region of interest, and detect peaks of PICs.
  For a single sample:

    library(KPIC)
    filename <- 'E:/LC-MS data/example/1-1_Seg1Ev1.mzXML'
    raw <- LoadData(filename)
    pics <- getPIC(raw, level=50000)
    # pics <- getPIC.kmeans(raw, level=50000)  ## use k-means 
    pics <- PICsplit(pics) # Optional
    pics <- getPeaks(pics)
    viewPICs(pics)

  For a set of samples:

    files <- 'E:/LC-MS data/example'
    PICS <- PICset(files, level=50000, export=F, par=T)
    # PICS <- PICset.kmeans(files, level=50000, export=F, par=T)  ## use k-means 
    PICS <- PICset.split(PICS)
    PICS <- PICset.getPeaks(PICS)
    viewPICs(PICS[[1]])


### Grouping and Alignment:
  Align PICs extracted in the last step across samples.

    groups_raw <- PICset.group(PICS, tolerance = c(0.1, 20))
    groups_align <- PICset.align(groups_raw, method='fftcc',move='loess')
    groups_align <- PICset.group(groups_align$picset, tolerance = c(0.1, 20))
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

    labels <- c(rep('A',6), rep('B',6)) # the class of each sample
    analyst.RF(labels, data$data.mat)

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
