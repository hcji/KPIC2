\name{PICset.fit}
\alias{PICset.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The \code{PICfit} function for a set of samples.
}
\description{
This function is used to apply \code{PICfit} method to the result of \code{PICresolve} founction.
}
\usage{
PICset.fit(picset, iter = 50, par = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{picset}{
The result object of \code{PICset.resolve} function.
}
  \item{iter}{
see \code{PICfit}.
}
  \item{par}{
Whether parallel method is used.
}
}

\value{
The processed \code{picset} object
}
