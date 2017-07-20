\name{PICset.split}
\alias{PICset.split}

\title{
The \code{PICsplit} function for a set of samples.
}
\description{
This function is used to apply \code{PICsplit} method to the result of \code{PICset} or \code{PICset.kmeans} founction.
}
\usage{
PICset.split(picset, par = FALSE)
}

\arguments{
  \item{picset}{
The result of \code{PICset} or \code{PICset.kmeans} founction.
}
  \item{par}{
Whether parallel method is used.
}
}

\value{
The processed \code{picset} object
}
