\name{KPICset}
\alias{KPICset}

\title{
Extracting the pure ion chromatograms for a set of samples
}
\description{
Return a KPICset object
}
\usage{
KPICset(files, roi_ppm=30, level = 500, alpha = 1, itol = 0.5,
                 min_snr = 3, min_ridge = 3, fst = 0.3, missp = 5, eval
                 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
path names of the files to read
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
A KPIC object
}

\author{
HC Ji
}