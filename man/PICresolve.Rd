\name{PICresolve}
\alias{PICresolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify whether overlapping peaks or not
}
\description{
With peak detection results, there are overlapping peaks in some mass trace. This function use the m/z values to identify whether they are real overlapping peaks or single peak with poor shape.
}
\usage{
PICresolve(pics, pval = 0.01)
}

\arguments{
  \item{pics}{
The result object of \code{split} function.
}
  \item{pval}{
The threshold value of m/z difference test.
}
}

\value{
a LIST of:
\item{scantime}{The retention time of each scan.}
\item{pics}{The extracted mass trace.}
\item{peaks}{The detected peak of each mass trace.}
}
