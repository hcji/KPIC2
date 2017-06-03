\name{KPICset}
\alias{KPICset}

\title{
Extracting the pure ion chromatograms for a set of samples
}
\description{
Return a KPICset object
}
\usage{
KPICset(files, roi_range = 0.1, level = 500, itol = c(-1, 1), min_snr = 3, peakwidth = c(5, 60),
        min_ridge =  3, fst = 0.3, missp = 5, cluster.ref = "square", eval = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
path names of the files to read
}
  \item{roi_range}{
numeric, the m/z range of ROI to be clursted
}
  \item{level}{
numeric, the minimum intensity of a feature, i.e. the minimum intensity of the landmark ion to be considered
}
  \item{itol}{
the relative tolerence between the intensity of ion signal and the estimated intensity, given as range (min,max). it could be Inf if estimated intensity is not used.
}
  \item{min_snr}{
numeric, the minimum SNR used for peak detection
}
  \item{peakwidth}{
Chromatographic peak width, given as range (min,max) in seconds
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
A KPIC object
}

\author{
HC Ji
}