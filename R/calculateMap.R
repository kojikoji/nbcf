##' Calculation of map estimates for change points
##'
##' This calculation is based on log probability for each observation intervals. see \code{\link{calculateLpst}}.
##' @title calculateMap
##' @param lpst Numeric matrix, log probability for each observation intervals
##' @param lambda Numeric, prior probability for change points in each transition
##' @return map.change.point, Numeric vector, MAP estimates of change points 
##' @author Yasuhiro Kojima
##' @import purrr
calculateMap <- function(lpst, lambda){
  ## calculate statics
  time.num <- ncol(lpst)
  end.stats <- list(lql = lpst[time.num, time.num], shatl = c(time.num))
  stats <- purrr::reduce(seq(time.num - 1, 1),
                             function(pre.stats, t){
                               lrl <- purrr::map(
                                               seq(t, time.num),
                                               function(s){
                                                 calculateLr(t, s, pre.stats$lql, lpst, lambda, time.num)
                                               }
                                             )
                               list(
                                 lql = c(max(unlist(lrl)), pre.stats$lql),
                                 shatl = c(which.max(unlist(lrl)) + (t - 1), pre.stats$shatl)
                               )
                             },
                             .init = end.stats
                             )
  ## trace back the shatl
  map.change.point <- traceGetShatl(stats$shatl, c(0))
}
##' Calculate statics for MAP estimate
##'
##' Calculate the log probability assuming  there is no change points between t ~ s, and change points are MAP for s+1 ~ time.num
##' @title calculateLr
##' @param t Integer, index of initial time points of unchange interval
##' @param s Integer, index of end time points of unchange interval
##' @param lql Numeric vector, log probability of interval t ~ time.num when change points are MAP
##' @param lpst Numeric matrix, log probability for each observation intervals
##' @param lambda Numeric, prior probability for change points in each transition
##' @param time.num Integer, index of end time points
##' @return lr Numeric, statics for MAP estimate
##' @author Yasuhiro Kojima
calculateLr <- function(t, s, lql, lpst, lambda, time.num){
  lpst[t, s] +  (s - t) * log(1 - lambda) + ifelse(s != time.num, log(lambda) + lql[s + 1 - t], 0)
}

##' Get trace of shatl
##'
##' This recursively get the elements of MAP change points, like \code{next.tau = shatl[tau + 1]}
##' @title traceGetShatl
##' @param shatl Numeric vector, next MAP change points of each index, the end end elements should be -1
##' @param taul Numeric vector, previous MAP change points
##' @return next.taul Numeric vector, one change points added for taul unless \code{shatl[taul[1] + 1] > 0}
##' @author 小嶋泰弘
traceGetShatl <- function(shatl, taul){
  if(taul[1] < max(shatl)){
    traceGetShatl(shatl, c(shatl[taul[1] + 1], taul))
  }else{
    rev(taul[taul != 0 & taul != max(shatl)])
  }
}
