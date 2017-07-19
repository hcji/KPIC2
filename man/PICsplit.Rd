\name{PICsplit}
\alias{PICsplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spliting multiple-peak trace into single ones
}
\description{
If there is more than one peak in a mass trace, and they are obviously separated, they can be split with this function.
}
\usage{
PICsplit(pics)
}

\arguments{
  \item{pics}{
The result of \code{getPIC} or \code{getPIC.kmeans} function.
}
}

\value{
a LIST of:
\item{scantime}{The retention time of each scan.}
\item{pics}{The extracted mass trace.}
\item{peaks}{The detected peak of each mass trace.}
}
