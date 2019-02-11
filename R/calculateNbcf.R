##' Function for detecting change points and changed variates
##'
##' This function is responsible for main calculation of this package. This calculate simulated change points and MAP estimates of them, and detect variates which significantly changed between each MAP estimated change points.
##' 
##' @title calculateNbcf
##' @param count.mat dgCmatrix, Count matrix. Each row represent each variate. Each column represent each observation
##' @param t.vec Numeric vector, Observation time
##' @param mean.bias.vec Numeric vector, Relative mean bias for each observation 
##' @param r Numeric, Size parameters of negative binomial distribution. 
##' @param lambda Numeric, Prior probability the change point occur each step
##' @param p.res Integer, Approximation resolution of parameter p
##' @param count.res Integer, Approximation resolution of count.
##' @param t.res Integer, Approximation resolution of t.
##' @return Nbcf, nbcf Instance of class \code{\link{Nbcf}}
##' @seealso [Nbcf]
##' @author Yasuhiro Kojima
##' 
##' @export

calculateNbcf <- function(count.mat, t.vec, mean.bias.vec,
                          alpha=1.0, beta=1.0, r=30, lambda=0.01,
                          p.res=1000, count.res=1000, t.res=100){
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec, r, alpha, beta, count.res, p.res)
  ## t.grids will seprate observation into t.res bins
  t.grids <- as.integer(seq(0, length(t.vec), length.out = t.res+1))
  ## count.mat is ordered based on t.vec 
  count.mat <- count.mat[, order(t.vec)]
  lpst <- calculateLpst(lnb.dict, count.mat, t.grids)
  lqt <- calculateLqt(lpst, lambda)
  indx.sim.change.point.list <- perfectSimulation(lqt, lpst, lambda)
  ## conver from grid index to t corresponding to end point of each grid
  sim.change.point.list <- purrr::map(
                                    indx.sim.change.point.list,
                                    ~ t.grids[.x + 1])
  map.change.point <- vector() ## estimateMap(lqt, lpst, lambda)
  change.variate <- list() ## detectVariate(map.change.point, count.mat, lnb.dict)
  nbcf <- new("Nbcf",
              count.mat = count.mat,
              t.vec = t.vec,
              mean.bias.vec = mean.bias.vec,
              params = list(
                alpha = alpha,
                beta = beta,
                r = r,
                lambda = lambda,
                p.res = p.res,
                count.res = count.res),
              lnb.dict = lnb.dict,
              lpst = lpst, lqt = lqt,
              sim.change.point.list = sim.change.point.list,
              map.change.point = map.change.point,
              change.variate = change.variate)
  return(nbcf)
}
