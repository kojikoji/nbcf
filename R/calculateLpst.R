##' Calculate the probability of count.mat[:, s:t]
##'
##' Based on probability values of grid points in \code{\link{lnb.dict}}, this calculate the probability values of specific time intervals of  \code{count.mat}
##' @title calculatePst
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
##' @title calculatePst
##' @param lnb.dict see \code{\link{lnb.dict}}
##' @param count.mat Integer matrix, count data. each row and column represent variate and observation
##' @param t.grids Numeric vector, grid index of t.vec
##' @return lpst Numeric matrix, log probability for each time interval
##' @author Yasuhiro Kojima
calculateLpstList <- function(lnb.dict, count.mat, t.grids){
  rcpp_calculate_lpst_list(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                           lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}


##' Calculate ratio of probability between divided and merged for local bin
##'
##' Based on probability for each intervals, comparing probability for fixed length interval which is divided at half point or not.
##' @title calculatePst
##' @param lpst Numeric matrix, log probability for each time interval
##' @return lhr.vec log ratio of compared probabilities
##' @author Yasuhiro Kojima
calculateLhrDivideTwo <- function(lpst, bin.prop=0.2){
  half.width <- as.integer((ncol(lpst)*bin.prop)/2)
  init <- half.width + 1
  end <- ncol(lpst) - (half.width - 1)
  lhr.vec <- unlist(purrr::map(seq(init, end),
                               ~ lpst[.x - half.width,(.x - 1)] + lpst[.x,.x + (half.width - 1)] - lpst[.x - half.width, .x + (half.width - 1)]
                               )
                    )
  return(lhr.vec)
}
