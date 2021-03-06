% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perfectSimulation.R
\name{decideSignificantChange}
\alias{decideSignificantChange}
\title{decideSignificantChange}
\usage{
decideSignificantChange(sim.change.point.list, sig.ct.thresh = 0.9)
}
\arguments{
\item{sim.change.point.list}{List, list of simulated change points}

\item{sig.ct.thresh}{Numeric, the proportion which detected point change appearance exceeds}
}
\value{
sig.ct
}
\description{
Detect change points which is repeatedly appears in perfectSimulation
}
\details{
This is instant version of detect ct
}
\author{
Yasuhiro Kojima
}
