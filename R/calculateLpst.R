##' Calculate the probability of count.mat[:, s:t]
##'
##' Based on probability values of grid points in \code{lnb.dict}, this calculate the probability values of specific time intervals of  \code{count.mat}
##' @title calculatePst
##' @param lnb.dict 
##' @param count.mat 
##' @param t.grids 
##' @return pst
##' @author 小嶋泰弘
calculateLpst <- function(lnb.dict, count.mat, t.grids){
  rcpp_calculate_lpst(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                      lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}

##' Calculate the probability of count.mat[:, s:t] variate-wise
##'
##' Based on probability values of grid points in \code{lnb.dict}, this calculate the probability values of specific time intervals of  \code{count.mat}. This returns the list for each variate
##' @title calculatePst
##' @param lnb.dict 
##' @param count.mat 
##' @param t.grids 
##' @return pst
##' @author 小嶋泰弘
calculateLpstList <- function(lnb.dict, count.mat, t.grids){
  rcpp_calculate_lpst_list(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                           lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
}
