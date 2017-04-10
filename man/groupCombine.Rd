\name{groupCombine}
\alias{groupCombine}

\title{
Combine tailed, isotopic or/and adduct features into the same group.
}
\description{
Combine tailed, isotopic or/and adduct features into the same group. 
}
\usage{
groupCombine(r.group,min_corr=0.9,type='tailed',window=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.group}{
the result of PICgroup function
}
  \item{min_corr}{
the minimum coefficient between peaks, which are regarded as isotopes or adducts and the base feature. It is recommended run min_corr = 0 for tailed features first time. and then run min_corr = 0.9 for isotopic features or adduct features. 
}
  \item{type}{
'tailed' for tailed features; 'isotope' for tailed features and isotopic features; or 'all'.
}
  \item{window}{
the width of RT window.
}
}
