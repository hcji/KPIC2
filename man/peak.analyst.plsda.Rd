\docType{methods}
\name{peak.analyst.plsda}
\alias{peak.analyst.plsda}
\title{peak.analyst.plsda}
\description{
PLS-DA classification, and cross-validation with the NIPALS algorithm
}

\usage{
peak.analyst.plsda(x, y, ncomp)
}

\arguments{
\item{x}{- n x p matrix, where n is the number of samples and p is the number of variables.}

\item{y}{- n x 1 matrix.}

\item{ncomp}{- Integer: number of components (predictive componenents in case of PLS.}
}
\value{
An object of class 'opls' (see ropls package) and VIPs
}


