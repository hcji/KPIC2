\name{getDataMatrix}
\alias{getDataMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate data matrix.
}
\description{
Extract all information from an KPIC object. Returns a peaklist.
}
\usage{
getDataMatrix(r.group, std = "maxo", frac=0.56)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.group}{
the result of grouping step
}
  \item{std}{
The standard for quantification, which could be 'maxo','peak_area'.
}
  \item{frac}{
the minimum fraction of samples which contain the same feature. 
}
}
