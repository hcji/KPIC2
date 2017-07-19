\name{LoadData}
\alias{LoadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load an LC-MS data file.
}
\description{
This function handles the task of reading a NetCDF/mzXML file containing LC-MS data.
}
\usage{
LoadData(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The path of LC-MS data file
}
}

\value{
A LIST of:
\item{mz}{The vector of m/z values.}
\item{scans}{The vector of scan indexes.}
\item{ints}{The vector of intensity values.}
\item{times}{The vector of unique time points.}
}

