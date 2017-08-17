\name{groupCombine}
\alias{groupCombine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine tailed, isotopic or/and adduct features into the same group.
}
\description{
Combine tailed, isotopic or/and adduct features into the same group.
}
\usage{
groupCombine(groups, min_corr = 0.9, type = "tailed", window = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{groups}{
The result of \code{PICset.group} function.
}
  \item{min_corr}{
the minimum coefficient between peaks, which are regarded as isotopes or adducts and the base feature.
}
  \item{type}{
'tailed' for tailed features; 'isotope' for tailed features and isotopic features; or 'all'.
}
  \item{window}{
the width of RT window.
}
}

\value{
a LIST of:
\item{peakmat}{The peakmat with group index.}
\item{picset}{The picset.}
\item{group.info}{The information of each group.}
}
