\name{getPIC.kmeans}
\alias{getPIC.kmeans}

\title{
Extract PICs from a LC-MS raw object based on optimal k-means clustering.
}
\description{
This method bases on the optimal k-means clustering of m/z values of data points in ROI. see reference for details.
}
\usage{
getPIC.kmeans(raw, level, mztol = 0.1, gap = 3, width = c(5, 60), alpha = 0.3, min_snr = 4, ...)
}

\arguments{
  \item{raw}{
Raw LC-MS data object obtained by \code{LoadData} function.
}
  \item{level}{
Mass traces are only retained if their maximum values are over \code{level}.
}
  \item{mztol}{
The m/z range of ROI.
}
  \item{gap}{
The number of gap points of a mass trace.
}
  \item{width}{
The range of a mass trace.
}
  \item{alpha}{
The parameter of forecasting.
}
  \item{min_snr}{
Minimum signal to noise ratio.
}
  \item{\dots}{
No use at present.
}
}

\value{
a LIST of:
\item{scantime}{The retention time of each scan.}
\item{pics}{The extracted mass trace.}
\item{peaks}{The detected peak of each mass trace.}
}
\references{
Ji, H., et al. "KPIC2: An Effective Framework for Mass Spectrometry-Based Metabolomics Using Pure Ion Chromatograms." Analytical Chemistry (2017).
}

\seealso{
\code{getPIC}
}
