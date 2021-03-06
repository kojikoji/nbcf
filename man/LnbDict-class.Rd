% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LnbDict.R
\docType{class}
\name{LnbDict-class}
\alias{LnbDict-class}
\title{Lnbdict}
\description{
Data stroge for log negative binomial probability
}
\details{
This class is a storage for precomputated log negative binomial probability. This also has the index conversion matrix and vector for parameter p and count.
}
\section{Slots}{

\describe{
\item{\code{values}}{Numeric Matrix, This storage log negative binomial probability. Each column and row represent each grid of count and p.}

\item{\code{lprior.values}}{Numeric vector, log prior probability values for each gridpoints}

\item{\code{count.dict}}{This is a sparse vector each elements have the row index of \code{value} cressponding the count represented by index of this vector}

\item{\code{p.dict}}{This is a dense matrix. Each element has a column index of \code{value} after correction against \code{\link[Nbcf]{mean.bias.vec}}. Each row represents index of p before correction. Each column represents observation index.}

\item{\code{p.width.vec}}{grid width of p}
}}

\author{
Yasuhiro Kojima
}
