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
PICgroup(xset,tolerance=c(0.01,10),weight=c(0.7,0.2,0.1),method='dbscan',frac=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
a KPIC object after alignment
}
  \item{tolerance}{
numeric, maximum allowed absolute m/z and RT difference
}
  \item{weight}{
numeric,  the assigned weight for m/z, RT and intensity difference at the moment of match score calculation between peaks.
}
  \item{method}{
'dbscan' or 'score'. 'dbscan' means group features with dbscan clustering method; 'score' means group features with calculated scores.
}
  \item{frac}{
minimum fraction of samples necessary in at least one of the sample groups.
}
}
