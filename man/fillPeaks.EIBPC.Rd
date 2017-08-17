\name{fillPeaks.EIBPC}
\alias{fillPeaks.EIBPC}

\title{
identify missing peaks
}
\description{
For each sample, identify missing peaks resulting from peak detection or other steps. The EIBPC is used to achieve this aim.
}
\usage{
fillPeaks.EIBPC(groups, extand_mz=20, extand_rt=5, min_snr=3, std='maxo')
}

\arguments{
  \item{groups}{
The result of \code{getDataMatrix} function
}
  \item{extand_mz}{
PPM of m/z tolerance of filled peaks.
}
  \item{extand_rt}{
Retention time tolerance of filled peaks.
}
  \item{min_snr}{
The minimum SNR of peaks, which may be smaller than that of getPIC function
}
  \item{std}{
The standard for quantification, only 'maxo' is supported now.
}
}
