\name{PICset.align}
\alias{PICset.align}

\title{
Align each group of PICs.
}
\description{
This function is used to calculated the shifts of PICs in each group, and correct the retention times of the peakmat and picset obtain by \code{PICset.group} function.
}
\usage{
PICset.align(groups, method = "fftcc", move = "direct", span = 1.5)
}

\arguments{
  \item{groups}{
The result of \code{PICset.group} function.
}
  \item{method}{
Which method is used to calculated the shift. can be \code{'match'} of \code{'fftcc'}. \code{'match'} means calculating the difference of the retention time of detected peak position. \code{'fftcc'} means use fft cross correlation method to maximize the similarity of peak shape.
}
  \item{move}{
Which method is used to move the original to new position. can be 'direct' or 'loess'. 'direct' means directly move each PIC based on the calculated shift. 'loess' means use a loess regression to the obtained shift and predict a new shift of each PIC, then move each PIC based on the new values.
}
  \item{span}{
The parameter which controls the degree of smoothing. Only used when the \code{move} is \code{'loess'}
}
}

\value{
a LIST of:
\item{peakmat}{The peakmat with refreshed \code{rt}.}
\item{picset}{The picset of refreshed \code{rt}.}
}
