##' Calculate Lqt, which represent probability of Yt:n given there is change point in t-1 
##'
##' These calculation are based on likelihood of each observation intervals
##' @title Calculate Lqt
##' @param lpst Numeric matrix, Log probability for each intervals
##' @param lambda Numeric, Prior probability for change point occurrence per step
##' @return lqt Numeric vector, Each value represents Yt:n probability
##' @author Yasuhiro Kojima
##' @export

calculateLqt <- function(lpst, lambda){
  t_len <- ncol(lpst)
  lqt <- vector()
  lqt[t_len] <- lpst[t_len, t_len]
  logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  calculateLqtVec <- function(t, s, pre.lqt){
    (s - t) * log(1 - lambda) + lpst[t,s] + ifelse(s != t_len, pre.lqt[s + 1 - t] + log(lambda), 0)
  }
  purrr::reduce(
           seq(t_len - 1, 1),
           function(pre.lqt, t){
             lqt_vec <- unlist(
               purrr::map(
                        t:t_len,
                        ~ calculateLqtVec(t, .x, pre.lqt)
                      )
             )
             c(logSumExp(lqt_vec), pre.lqt)
           },
           .init = lpst[t_len, t_len]
         )
}
