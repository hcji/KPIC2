\name{PICfit}
\alias{PICfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit peaks with Gaussian curve.
}
\description{
Fit all of the detected peaks with Gaussian cure, including the overlapping peaks.
}
\usage{
PICfit(pics, iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pics}{
The result object of \code{PICresolve} or \code{PICsplit} function.
}
  \item{iter}{
Maximum iterations.
}
}

\value{
a LIST of:
\item{scantime}{The retention time of each scan.}
\item{pics}{The fitted curve.}
\item{peaks}{The detected peaks.}
}
