\name{getPIC}
\alias{getPIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting the pure ion chromatograms via Optimal K-means Clustering
}
\description{
Return a list which is comprised of PICs and the main information of them.
}
\usage{
getPIC(filename, roi_ppm = 30, level = 500, alpha = 1, itol = 0.5, min_snr = 3, min_ridge = 3, fst = 0.3, missp = 5)
}
\arguments{
  \item{filename}{
char, path of the LC-MS data file to read, can be.mzXML, .CDF, .mzdata, .mzml, etc
}
  \item{roi_ppm}{
numeric, the m/z range of ROI to be clursted
}
  \item{level}{
numeric, the minimum intensity of a feature, i.e. the minimum intensity of the landmark ion to be considered
}

  \item{itol}{
numeric, the relative tolerence between the intensity of ion signal and the estimated intensity. it could be Inf if estimated intensity is not used.
}
  \item{min_snr}{
numeric, the minimum SNR used for peak detection
}
  \item{min_ridge}{
numeric, the minimum scale of peaks
}
  \item{fst}{
numeric, the parameter for intensity forcasting via exponential smoothing.
}
  \item{missp}{
numeric, number of missing point allowed in a PIC.
}
  \item{eval}{
logistic, whether to evaluate the peak shape or not.
}
}

\value{
a KPIC list of:
  \item{PICs}{
a list of pure ion chromatograms
}
  \item{Info}{
  the information of each PIC:
  
  \code{mz} mean m/z of each PIC
    
  \code{mzmin} minimum m/z
    
  \code{mzmax} maximum m/z
    
  \code{rt} position of peak
    
  \code{rtmin} left boundary of peak
    
  \code{rtmax} right boundary of peak
    
  \code{maxo} maximum intensity
    
  \code{peak_area} integrated peak area
    
  \code{snr} SNR of the peak
    
    \code{scale} the fitted scale of wavelet
    
    \code{wt_coef} the fitted wavelet coefficient
}
}

\author{
HC Ji
}
