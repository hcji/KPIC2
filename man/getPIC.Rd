\name{getPIC}
\alias{getPIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract PICs from a LC-MS raw object based on m/z difference.
}
\description{
This method bases on the extension of mass trace depending on the m/z difference. The tolerence is described via \code{mean} and \code{variance}.
}
\usage{
getPIC(raw, level, mztol = 0.1, gap = 3, width = 5, min_snr = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
Raw LC-MS data object obtained by \code{LoadData} function.
}
  \item{level}{
Mass traces are only retained if their maximum values are over \code{level}.
}
  \item{mztol}{
The initial m/z tolerence.
}
  \item{gap}{
The number of gap points of a mass trace.
}
  \item{width}{
The minimum length of a mass trace.
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

\seealso{
\code{getPIC.kmeans}
}
