\name{fillPeaks.peakfinder}
\alias{fillPeaks.peakfinder}

\title{
identify missing peaks
}
\description{
For each sample, fill missing peaks through find peaks form KPICset. 
}
\usage{
fillPeaks.peakfinder(r.DataMatrix,tolerance=c(0.1,15),weight=c(0.7,0.2,0.1),std='maxo')
}

\arguments{
  \item{r.DataMatrix}{
The result of \code{getDataMatrix} function
}
  \item{tolerance}{
The tolerance of filled peaks: c(m/z difference, rt difference)
}
  \item{weight}{
The weight of m/z difference, rt difference and intensities difference
}
  \item{std}{
The standard for quantification, which could be 'maxo','peak_area'.
}
}
