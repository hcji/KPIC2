\name{PICgroup}
\alias{PICgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
group the features
}
\description{
group the features across samples
}
\usage{
PICgroup(xset, beta = 0.8, group.width = 15, ppm = 80, min_CoP = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
a KPIC object after alignment
}
  \item{beta}{
numeric, a factor of group width, which is calculated by beta*FWHM of base feature.
}
  \item{group.width}{
numeric, the width of each group. if NULL, it will be calculated by beta; if given, the beta is invalid.
}
  \item{ppm}{
the tolerance of m/z of features belonging to the same group
}
  \item{min_CoP}{
the minimum coefficient between peaks and the base feature, which are regarded as isotopes or adducts.
}
}
