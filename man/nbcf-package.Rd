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
data(sim.data)
nbcf.obj <- calculateNbcf(count.mat=sim.data$count.mat, t.vec=sim.data$t.vec, alpha=1.0, beta=1.0, r=30, lambda= 0.01)
hist(nbcf.obj@change.point)
nbcf.obj@change.variate
}
