##' Calculate Lqt, which represent probability of Yt:n given there is change point in t-1 
##'
##' These calculation are based on likelihood of each observation intervals
##' @title Calculate Lqt
##' @param lpst Numeric matrix, Log probability for each intervals
##' @param lambda Numeric, Prior probability for change point occurrence per step
##' @return lqt Numeric vector, Each value represents Yt:n probability
##' @author Yasuhiro Kojima
calculateLqt <- function(lpst, lambda){
  t_len <- ncol(lpst)
  lqt <- vector()
  lqt[t_len] <- lpst[t_len, t_len]
  logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  for(t in seq(t_len - 1, 1)){
    lqt_vec <- unlist(purrr::map(
                               t:(t_len-1),
                               ~ lqt[.x + 1] + log(lambda) + (.x - t) * log(1 - lambda) + lpst[t,.x]))
    no_change_lqt <- (t_len - t) * log(1 - lambda) + lpst[t, t_len]
    lqt_vec <- c(lqt_vec, no_change_lqt)
    lqt[t] <- logSumExp(lqt_vec)
  }
  return(lqt)
}
