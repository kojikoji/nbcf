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
##' @param map.fix.num Integer, the number of map change points. If it is \code{NULL}, the point number is automaticaly estimated
##' @return Nbcf, nbcf Instance of class \code{\link{Nbcf}}
##' @seealso [Nbcf]
##' @author Yasuhiro Kojima
##'
##' @import purrr
##' @export

calculateNbcf <- function(count.mat, t.vec, mean.bias.vec,
                          alpha=1.0, beta=1.0, r=30, lambda=0.01,
                          p.res=1000, count.res=1000, t.res=100, map.fix.num = NULL){
  ## count.mat and mean.bias.vec  are ordered based on t.vec 
  count.mat <- count.mat[, order(t.vec)]
  mean.bias.vec <- mean.bias.vec[order(t.vec)]
  t.vec <- t.vec[order(t.vec)]
  ## rows of count.mat will be named as seq(nrow(count.mat))
  if(length(rownames(count.mat)) == 0){
    rownames(count.mat) <- seq(nrow(count.mat))
  }
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec, r, alpha, beta, count.res, p.res)
  ## t.grids will seprate observation into t.res bins
  t.grids <- as.integer(seq(0, length(t.vec), length.out = t.res+1))
  lpst.list <- calculateLpstList(lnb.dict, count.mat, t.grids)
  names(lpst.list) <- rownames(count.mat)
  lqt.list <- purrr::map(lpst.list, ~ calculateLqt(.x, lambda))
  bf.list <- purrr::map2(lqt.list, lpst.list, ~ .x[1] - .y[1, ncol(.y)])
  lpst <- purrr::reduce(lpst.list, ~ .x + .y)
  used.vars <- rownames(count.mat)
  lqt <- calculateLqt(lpst, lambda)
  sim.change.point.list <- perfectSimulation(lqt, lpst, lambda)
  if(is.numeric(map.fix.num)){
    map.change.point <- calculateMapFix(lpst, lambda, map.fix.num)
  }else{
    map.change.point <- calculateMap(lpst, lambda)
  }
  change.variate.df <- findChangeVariate(lpst.list, used.vars, map.change.point, lambda = lambda)
  ## conver from grid index to t corresponding to end point of each grid
  convertTidx2T <- function(t.idx.vec){
    unlist(purrr::map(t.idx.vec, ~ t.vec[t.grids[.x + 1]]))
  }
  outer.sim.change.point.list <- purrr::map(
                                    sim.change.point.list,
                                    ~ t.vec[t.grids[.x + 1]])
  outer.map.change.point <- convertTidx2T(map.change.point)
  change.variate.df$change.point <- convertTidx2T(change.variate.df$change.point)
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
              lpst = lpst,
              lqt = lqt,
              lpst.list = lpst.list,
              lqt.list = lqt.list,
              bf.list = bf.list,
              used.vars = used.vars,
              sim.change.point.list = outer.sim.change.point.list,
              map.change.point = outer.map.change.point,
              change.variate.df = change.variate.df)
  return(nbcf)
}
