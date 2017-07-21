\name{PICset.group}
\alias{PICset.group}

\title{
group the features
}
\description{
This function is used to group the features across samples.
}
\usage{
PICset.group(picset, tolerance = c(0.01, 10), weight = c(0.8, 0.2), method = "score", frac = 0.5)
}

\arguments{
  \item{picset}{
The result of \code{PICset.getPeaks} function.
}
  \item{tolerance}{
Maximum allowed absolute m/z and RT difference
}
  \item{weight}{
The assigned weight for m/z and RT difference at the moment of match score calculation between peaks.
}
  \item{method}{
Which method is used. can be 'score' or 'dbscan'. 'dbscan' means group features with dbscan clustering method; 'score' means group features with calculated scores.
}
  \item{frac}{
Minimum fraction of samples necessary in at least one of the sample groups.
}
}

\value{
a LIST of:
\item{peakmat}{The final peakmat of all sample with group id.}
\item{picset}{The input picset}
}
