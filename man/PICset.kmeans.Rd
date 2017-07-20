\name{PICset.kmeans}
\alias{PICset.kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process a set of sample with \code{getPIC.kmeans} method.
}
\description{
This function is used to process a dataset produced by LC-MS.
}
\usage{
PICset.kmeans(files, level, mztol = 0.1, gap = 3, width = c(5, 60), min_snr = 4, alpha = 0.3, equal = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
The path of the LC-MS files folder.
}
  \item{level}{
see \code{getPIC.kmeans}
}
  \item{mztol}{
see \code{getPIC.kmeans}
}
  \item{gap}{
see \code{getPIC.kmeans}
}
  \item{width}{
see \code{getPIC.kmeans}
}
  \item{min_snr}{
see \code{getPIC.kmeans}
}
  \item{alpha}{
see \code{getPIC.kmeans}
}
  \item{equal}{
Whether the retention times of samples are equaled or not. Equalization is need for alignment procedure.
}
  \item{\dots}{
see \code{getPIC.kmeans}
}
}

\value{
a LIST of PICs, each element is the result of \code{getPIC.kmeans} function.
}
