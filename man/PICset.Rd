\name{PICset}
\alias{PICset}

\title{
Process a set of sample with \code{getPIC} method.
}
\description{
This function is used to process a dataset produced by LC-MS.
}
\usage{
PICset(files, level, mztol = 0.1, gap = 3, width = 5, min_snr = 4, equal = TRUE, export=FALSE, par=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
The path of the LC-MS files folder.
}
  \item{level}{
see \code{getPIC}
}
  \item{mztol}{
see \code{getPIC}
}
  \item{gap}{
see \code{getPIC}
}
  \item{width}{
see \code{getPIC}
}
  \item{min_snr}{
see \code{getPIC}
}
  \item{equal}{
Whether the retention times of samples are equaled or not. Equalization is need for alignment procedure.
}
  \item{export}{
Whether to export PICs of each sample as single files
}
  \item{par}{
Whether to use multi-core calculation
}
  \item{\dots}{
see \code{getPIC}
}
}

\value{
a LIST of PICs, each element is the result of \code{getPIC} function.
}
