\name{getPeaks}
\alias{getPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get peaks of the detected PICs.
}
\description{
Get the information peaks of the detected PICs, including m/z, retention time, snr, scale, height and peak area, etc. Note, only the information of highest peak of a PIC will be included.
}
\usage{
getPeaks(pics)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pics}{
The result object of \code{getPIC}, \code{getPIC.kmeans}, \code{PICsplit}, \code{PICresolve} or \code{PICfit} founction.
}
}
\value{
\item{scantime}{The retention time of each scan.}
\item{pics}{The extracted mass trace.}
\item{peaks}{The detected peak of each mass trace.}
\item{peakInfo}{The information of the peaks.}
}
