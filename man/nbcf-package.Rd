\name{nbcf-package}
\alias{nbcf-package}
\alias{nbcf}
\docType{package}
\title{
  Change point detection for time course multi negative binomial observation
}
\description{
  This is a package for calculating change point for time course multi negative binomial observation. Furthermore, this package aims to detect the variate whose distribution is significantly changed across each change point. Coming soon.
}
\details{
  Comming soon
}
\author{
Yasuhiro Kojima, tiisaishima@gmail.com

Maintainer: Yasuhiro Kojima <tiisaishima@gmail.com>
}
\references{
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}

\examples{
sim.data <- data(sim.data)
nbcf.obj <- loadData(sim.data$count.mat, sim.data$t)
nbcf.obj <- setupNegativeBinomDict(nbcf.obj)
nbcf.obj <- calculateQt(nbcf.obj, lambda= 0.01, alpha=1.0, beta=1.0)
nbcf.obj <- perfectSimulation(nbcf.obj)
nbcf.obj <- detectVariate(nbcf.obj)
hist(nbcf.obj@change.point)
nbcf.obj@change.variate
}
