##' Calculate the probability of count.mat[:, s:t]
##'
##' Based on probability values of grid points in \code{\link{lnb.dict}}, this calculate the probability values of specific time intervals of  \code{count.mat}
##' @title calculateLpst
##' @param lnb.dict see \code{\link{lnb.dict}}
##' @param count.mat Integer matrix, count data. each row and column represent variate and observation
##' @param t.grids Numeric vector, grid index of t.vec
##' @return lpst Numeric matrix, log probability for each time interval
##' @author Yasuhiro Kojima
calculateLpst <- function(lnb.dict, count.mat, t.grids){
  rcpp_calculate_lpst(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                      lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}

##' Calculate the probability of count.mat[:, s:t] variate-wise
##'
##' Based on probability values of grid points in \code{lnb.dict}, this calculate the probability values of specific time intervals of  \code{count.mat}. This returns the list for each variate
##' @title calculateLpstList
##' @param lnb.dict see \code{\link{lnb.dict}}
##' @param count.mat Integer matrix, count data. each row and column represent variate and observation
##' @param t.grids Numeric vector, grid index of t.vec
##' @return lpst Numeric matrix, log probability for each time interval
##' @author Yasuhiro Kojima
calculateLpstList <- function(lnb.dict, count.mat, t.grids){
  rcpp_calculate_lpst_list(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                           lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}


##' Calculate the probability of count from  time poinst \code{s} to \code{t} for one variate
##'
##' Based on probability values of grid points in \code{lnb.dict}, this calculate the probability values of specific time intervals of  \code{count.vec}. This returns a matrix whose row and column represent initial and end points of each interval
##' @title calculateLpstOneVariate
##' @param lnb.dict see \code{\link{lnb.dict}}
##' @param count.vec Integer vector, count data. each elements represent each observation
##' @param t.grids Numeric vector, grid index of t.vec
##' @return lpst Numeric matrix, log probability for each time interval
##' @author Yasuhiro Kojima
calculateOneVariateLpst <- function(lnb.dict, count.vec, t.grids){
  rcpp_calculate_lpst_g(count.vec, lnb.dict@values, lnb.dict@lprior.values,
                        lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}


##' Calculate ratio of probability between divided and merged for local bin
##'
##' Based on probability for each intervals, comparing probability for fixed length interval which is divided at half point or not.
##' @title calculateLhrDivideTwo
##' @param lpst Numeric matrix, log probability for each time interval
##' @return lhr.vec log ratio of compared probabilities
##' @author Yasuhiro Kojima
##' @export

calculateLhrDivideTwo <- function(lpst, bin.prop=0.2){
  half.width <- as.integer((ncol(lpst)*bin.prop)/2)
  init <- half.width + 1
  end <- ncol(lpst) - (half.width - 1)
  t.vec <- seq(init, end)
  names(t.vec) <- as.character(t.vec)
  lhr.vec <- unlist(purrr::map(t.vec,
                               ~ lpst[.x - half.width,(.x - 1)] + lpst[.x, .x + (half.width - 1)] - lpst[.x - half.width, .x + (half.width - 1)]
                               )
                    )
  return(lhr.vec)
}


##' Calculate variance of log normalized count for each intervals and each variate
##'
##' Calculate variance of log normalized count for each intervals and each variate. This is used when \code{method == "variance"} in calculateNbcf.
##' @title calculateLpstVar
##' @param count.mat Integer matrix, count data. each row and column represent variate and observation
##' @param mean.bias.vec Numeric vector, bias for each observation
##' @param t.grids Numeric vector, grid index of t.vec
##' @return lpst.list List, variance of log normalized counts at each interval and each variate
##' @author Yasuhiro Kojima
##' @import purrr
calculateVarLpstList <- function(count.mat, mean.bias.vec, t.grids){
    ## functions set  
    normalizeLogNormalize <- function(count.vec){
      log(count.vec/mean.bias.vec + 1) %>%
	{length(.) * . / sum(.)}
    }
    calculateElement <- function(idx, norm.vec){
      begin.t <- t.grids[idx]
      end.t <- t.grids[idx + 1] - 1
      part.vec <- norm.vec[begin.t:end.t]
      N <- end.t - begin.t + 1
      mu <- sum(part.vec)/N
      V <- sum((part.vec - mu)**2)
      list(
	N = N,
	mu = mu,
	V = V
      )
    }
    calculateN <- function(pre.stats, current.stats){
      c(pre.stats$N, 0) + current.stats$N
    }
    calculateMu <- function(pre.stats, current.stats, N){
      (c(pre.stats$N * pre.stats$mu, 0) + current.stats$N * current.stats$mu)/N
    }
    calculateV <- function(pre.stats, current.stats, N, mu){
      c(pre.stats$V + pre.stats$N * (pre.stats$mu - mu[1:length(pre.stats$mu)])**2, 0) + current.stats$V + current.stats$N * (current.stats$mu - mu)**2
    }
    calculateColStats <- function(pre.stats, current.stats){
      N <- calculateN(pre.stats, current.stats)
      mu <- calculateMu(pre.stats, current.stats, N)
      V <- calculateV(pre.stats, current.stats, N, mu)
      list(
	N = N,
	mu = mu,
	V = V
      )
    }
    calculateVar <- function(var){
      norm.vec <- normalizeLogNormalize(count.mat[var,])
      element.stats.list <- purrr::map(seq(end.idx), calculateElement, norm.vec)
      col.stats.list <- purrr::accumulate(element.stats.list, calculateColStats)
      base.mat <- purrr::map(col.stats.list, ~ c(-.x$V, rep(0, end.idx - length(.x$V)))) %>%
	{do.call(cbind, .)}
      ## mod.mat <- purrr::map(col.stats.list, ~ c(.x$N/(.x$N - 1), rep(0, end.idx - length(.x$N)))) %>% {do.call(cbind, .)}
      base.mat ##  * mod.mat
    }
    ## parameter set
    end.idx <- length(t.grids) - 1
    lpst.list <- purrr::map(rownames(count.mat), calculateVar)
  return(lpst.list)
}
