\name{PIAlign}
\alias{PIAlign}
\title{Correct retention time from different samples}
\description{To correct differences between retention times between different samples}
\usage{PIAlign(xset, mzSegSize = 0.5, shift = 10, lambda = 1.5, ref = NA)}

\arguments{
  \item{xset}{a KPICSet object.}
  \item{mzSegSize}{width of m/z range included in a MSC}
  \item{shift}{maximum shift of a signal}
  \item{lambda}{the coefficient of smoother}
  \item{ref}{the index of the reference sample, if missing, it will be optimized automatically}
}

\value{
the KPICSet object.
}

\author{
HC Ji
}