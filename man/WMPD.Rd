\name{WMPD}
\alias{WMPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reslove overlapped peak based on mass spectrometry.
}
\description{
Reslove overlapped peak based on mass spectrometry.
}
\usage{
WMPD(pic, min_snr, level, pval, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pic}{
Extracted ion trace.
}
  \item{min_snr}{
Minimum signal to noise ratio.
}
  \item{level}{
Peaks are only retained if their maximum values are over \code{level}.
}
  \item{pval}{
The p-value threshold of different peaks.
}
  \item{iter}{
Number of iteration
}
}

\value{
The result of peak detection.
}
