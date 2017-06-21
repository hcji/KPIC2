\name{fillPeaks.EIBPC}
\alias{fillPeaks.EIBPC}

\title{
identify missing peaks
}
\description{
For each sample, identify missing peaks resulting from peak detection or other steps. The EIBPC is used to achieve this aim.
}
\usage{
fillPeaks(r.DataMatrix,tolerance=c(0.1,15),min_snr=3,min_ridge=2,std='maxo',peak='highest')
}

\arguments{
  \item{r.DataMatrix}{
The result of \code{getDataMatrix} function
}
  \item{xset}{
A KPIC object
}
  \item{tolerance}{
The tolerance of filled peaks: c(m/z difference, rt difference)
}
  \item{min_snr}{
The minimum SNR of filled peak, which may be smaller than that of getPIC function
}
  \item{min_ridge}{
The minimum scale of peaks, which may be smaller than that of getPIC function
}
  \item{std}{
The standard for quantification, only 'maxo' is supported now.
}
  \item{peak}{
The standard for miss peak searching of EIBPC, which could be 'nearest','highest'.
}
}
