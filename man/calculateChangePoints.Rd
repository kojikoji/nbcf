% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateChangePoints.R
\name{calculateChangePoints}
\alias{calculateChangePoints}
\title{calculateChangePoints}
\usage{
calculateChangePoints(count.mat, t.vec, mean.bias.vec, alpha = 1,
  beta = 1, r = 30, lambda = 0.01, p.res = 1000,
  count.res = 1000, t.res = 100, sim.iter = 100, method = "map",
  map.num = 1)
}
\arguments{
\item{count.mat}{dgCmatrix, Count matrix. Each row represent each variate. Each column represent each observation}

\item{t.vec}{Numeric vector, Observation time}

\item{mean.bias.vec}{Numeric vector, Relative mean bias for each observation}

\item{r}{Numeric, Size parameters of negative binomial distribution.}

\item{lambda}{Numeric, Prior probability the change point occur each step}

\item{p.res}{Integer, Approximation resolution of parameter p}

\item{count.res}{Integer, Approximation resolution of count.}

\item{t.res}{Integer, Approximation resolution of t.}

\item{sim.iter}{Integer, Iteration of perfect simulation. See \code{\link{perfectSimulation}}}

\item{method}{Specify the method to calculate change points. \code{"map"} or \code{"sim"}}

\item{map.num}{Number of change points. This works only if you specify method as \code{"map"}}
}
\value{
change.point.df, Data frame, This contains simulated change points for each variates. This is composed of columns of changed variates, change point coordinates.
}
\description{
Function for detecting change points and changed variates
}
\details{
This function is responsible for main calculation of this package. This calculate simulated change points for each variates, and provide them a data frame
}
\author{
Yasuhiro Kojima
}
