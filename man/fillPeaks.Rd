\name{fillPeaks}
\alias{fillPeaks}

\title{
identify missing peaks
}
\description{
For each sample, identify missing peaks resulting from peak detection or other steps. The EIBPC is used to achieve this aim.
}
\usage{
fillPeaks(r.DataMatrix,expand_rt=1.5,expand_mz=150,min_snr=3,min_ridge=2,std='maxo')
}

\arguments{
  \item{r.DataMatrix}{
The result of \code{getDataMatrix} function
}
  \item{xset}{
A KPIC object
}
  \item{expand_rt}{
Expansion factor for the rentention time range, which is the times of peak width.
}
  \item{expand_mz}{
Expansion factor for the m/z range, which is ppm of m/z.
}
  \item{min_snr}{
The minimum SNR of filled peak, which may be smaller than that of getPIC function
}
  \item{min_ridge}{
The minimum scale of peaks, which may be smaller than that of getPIC function
}
  \item{std}{
The standard for quantification, which could be 'maxo','peak_area'.
}
}
