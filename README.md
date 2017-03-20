# KPIC2
  An Effective Framework for Mass Spectrometry-Based Metabolomics Using Pure Ion Chromatograms

## Installation  

####Install Depends: 
		install.packages(c("devtools","Rcpp","Ckmeans.1d.dp","foreach","doParallel","iterators","randomForest", "forecast", "pracma", "caret", "ellipse", "splancs", "ggplot2", "gtools"))
		source("https://bioconductor.org/biocLite.R")
		biocLite("mzR","ropls")
		install "kopls" package from kopls_1.1.2.zip, which is embodied in KPIC2

####Install KPIC:  

		library(devtools);  
		install_github("hcji/KPIC2")

## Usage 
The user guide as well as example scripts and datasets are integrated in the package.
  
## Contact
  ji.hongchao@foxmail.com