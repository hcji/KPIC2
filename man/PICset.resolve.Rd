\name{PICset.resolve}
\alias{PICset.resolve}

\title{
The \code{PICresolve} function for a set of samples.
}
\description{
This function is used to apply \code{PICresolve} method to the result of \code{PICsplit} founction.
}
\usage{
PICset.resolve(picset, pval = 0.01, par = FALSE)
}

\arguments{
  \item{picset}{
The result object of \code{PICset.split} function.
}
  \item{pval}{
see \code{PICsplit}.
}
  \item{par}{
Whether parallel method is used.
}
}

\value{
The processed \code{picset} object
}
