##' Function for detecting change points and changed variates
##'
##' This function is responsible for main calculation of this package. This calculate simulated change points for each variates, and provide them a data frame
##' 
##' @title calculateChangePoints
##' @param count.mat dgCmatrix, Count matrix. Each row represent each variate. Each column represent each observation
##' @param t.vec Numeric vector, Observation time
##' @param mean.bias.vec Numeric vector, Relative mean bias for each observation 
##' @param r Numeric, Size parameters of negative binomial distribution. 
##' @param lambda Numeric, Prior probability the change point occur each step
##' @param p.res Integer, Approximation resolution of parameter p
##' @param count.res Integer, Approximation resolution of count.
##' @param t.res Integer, Approximation resolution of t.
##' @param sim.iter Integer, Iteration of perfect simulation. See \code{\link{perfectSimulation}}
##' @return change.point.df, Data frame, This contains simulated change points for each variates. This is composed of columns of changed variates, change point coordinates.
##' @author Yasuhiro Kojima
##'
##' @import tibble
##' @import dplyr
##' @export

calculateChangePoints <- function(count.mat, t.vec, mean.bias.vec,
                          alpha=1.0, beta=1.0, r=30, lambda=0.01,
                          p.res=1000, count.res=1000, t.res=100, sim.iter=100, min.lhr=10){
  ## count.mat and mean.bias.vec  are ordered based on t.vec 
  count.mat <- count.mat[, order(t.vec)]
  mean.bias.vec <- mean.bias.vec[order(t.vec)]
  t.vec <- t.vec[order(t.vec)]
  ## rows of count.mat will be named as seq(nrow(count.mat))
  if(length(rownames(count.mat)) == 0){
    rownames(count.mat) <- seq(nrow(count.mat))
  }
  ## calculate negative binomial probabilities for discretized p and count and
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec, r, alpha, beta, count.res, p.res)
  ## t.grids will seprate observation into t.res bins
  t.grids <- as.integer(seq(0, length(t.vec), length.out = t.res+1))
  ## variatewise change points estimation
  change.point.df <- purrr::map_dfr(
    rownames(count.mat),
    function(var){
      lpst <- calculateOneVariateLpst(lnb.dict, count.mat[var,], t.grids)
      lqt <- calculateLqt(lpst, lambda)
      bayes.factor <- lqt[1] - lpst[1, ncol(lpst)]
      map.change.point <- calculateMap(lpst, lambda)
      ## put NA for no change point variate
      if(length(map.change.point) == 0) map.change.point <- NA
      tibble(var = var, change.point = map.change.point, bayes.factor = bayes.factor)
    })
  ##:ess-bp-start::browser@nil:##
  change.point.df <- change.point.df %>%
    ## conver from grid index to t corresponding to end point of each grid
    dplyr::mutate(change.point = t.vec[t.grids[change.point + 1]]) %>%
    ## calculate FDR based on bayes factor
    dplyr::arrange(dplyr::desc(bayes.factor)) %>%
    dplyr::mutate(FDR = cumsum(1/(1 + exp(bayes.factor)))/dplyr::row_number(dplyr::desc(bayes.factor)))
  return(change.point.df)
}
