##' Function for detecting change points and changed variates
##'
##' This function is responsible for main calculation of this package. This calculate simulated change points and MAP estimates of them, and detect variates which significantly changed between each MAP estimated change points.
##' 
##' @title calculateNbcf
##' @param count.mat Count matrix. Each row represent each variate. Each column represent each observation
##' @param t.vec Observation time
##' @param mean.bias.vec Relative mean bias for each observation 
##' @param alpha Parameter fo beta distribution, which is prior for p of negative binomial distribution
##' @param beta Parameter fo beta distribution, which is prior for p of negative binomial distribution
##' @param r Size parameters of negative binomial distribution. 
##' @param lambda Prior probability the change point occur each step
##' @param p.res Approximation resolution of parameter p
##' @param count.res Approximation resolution of count.
##' @param t.res Approximation resolution of t.
##' @return nbcf Instance of class \code{\link{Nbcf}}
##' @seealso [Nbcf]
##' @author Yasuhiro Kojima

calculateNbcf <- function(count.mat, t.vec, mean.bias.vec,
                          alpha=1.0, beta=1.0, r=30, lambda=0.01,
                          p.res=1000, count.res=1000, t.res=100){
  lnb.dict <- setupLnbDict(count.mat, mean.bias.vec, p.res, count.res)
  Pst <- calculatePst(lnb.dict.mat, count.mat, alpha, beta, r, t.res)
  Qt <- calculateQt(Pst, lambda)
  sim.change.point <- perfectSimulation(Qt, Pst, lambda)
  map.change.point <- estimateMap(Qt, Pst, lambda)
  change.variate <- detectVariate(map.change.point, count.mat, lnb.dict)
  nbcf <- new("Nbcf", count.mat=count.mat, t.vec=t.vec, mean.bias.vec=mean.bias.vec,
                  params=list(alpha=alpha, beta=beta, r=r, lambda=lambda,
                              p.res=p.res, count.res=count.res),
                  lnb.dict=lnb.dict, Pst=Pst, Qt=Qt,
                  sim.change.point=sim.change.point, map.change.point=map.change.poin,
                  change.variate=change.variate)
  return(nbcf)
}
