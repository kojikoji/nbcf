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
      {./sum(.)}
  }
  calculateElement <- function(init.idx, row.idx, norm.vec){
    norm.vec[t.grids[row.idx]:t.grids[init.idx + 1]] %>%
      {-sum((. - mean(.))**2)}
  }
  calculateRow <- function(row.idx, norm.vec){
    c(rep(0, row.idx-1), unlist(purrr::map(row.idx:end.idx, calculateElement, row.idx, norm.vec)))
  }
  calculateVar <- function(var){
    norm.vec <- normalizeLogNormalize(count.mat[var,]) * 10000
    purrr::map(seq(end.idx), calculateRow, norm.vec) %>%
      {do.call(rbind, .)} 
  }
  ## parameter set
  end.idx <- length(t.grids) - 1
  lpst.list <- purrr::map(rownames(count.mat), calculateVar)
  return(lpst.list)
}
