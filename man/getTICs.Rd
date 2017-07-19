\name{getTICs}
\alias{getTICs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get TICs of LC-MS data.
}
\description{
Get TICs of LC-MS data.
}
\usage{
getTICs(files, method = "BPC")
}

\arguments{
  \item{files}{
The path of LC-MS files.
}
  \item{method}{
\code{TIC} or \code{BPC}
}
}

\value{
a LIST, use
\item{rt}{retention time of each scan.}
\item{tics}{obtained tics.}
}
