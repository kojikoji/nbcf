##' Find change variates at each change points
##'
##' This function calculate the relation ships between change points and variates which have one change points among previous and following change points for each change point. Furthermore, we calculate the entropy of the change point and following change point.
##' @title findChangeVariate
##' @param lpst.list list of Numeric matrix, list of log probability for each intervals for each variate
##' @param used.vars vector, names of variates used in change point finding
##' @param ct.vec Numeric vector, change point vector
##' @return change.var.df tibble, contains variate names, corresponding change point, entropy and change number among corresponding intervals
##' @author Yasuhiro Kojima
findChangeVariate <- function(lpst.list, used.vars, ct.vec, lambda = 0.5){
  ## add ct vec to begin and end of time
  extended.ct.vec <- c(0, ct.vec, ncol(lpst.list[[1]]))
  df <- crossing(var = used.vars, ct.idx = seq(length(ct.vec))) %>%
      dplyr::mutate(
             ct = purrr::map(ct.idx, ~ ct.vec[.x]),
             lpst = purrr::map2(
                             var, ct.idx,
                             ~ lpst.list[[.x]][(extended.ct.vec[.y - 1 + 1] + 1):extended.ct.vec[.y + 1 + 1],
                                               (extended.ct.vec[.y - 1 + 1] + 1):extended.ct.vec[.y + 1 + 1]]
                             ),
             lqt = purrr::map(
                            lpst,
                            ~ calculateLqt(.x, lambda)
                          ),
             sim.ct.list = purrr::map2(
                                    lqt, lpst,
                                    ~ perfectSimulation(.x, .y, lambda)
                                  ),
             median.change.num = median(
               unlist(
                 purrr::map(
                          sim.ct.list,
                          ~ length(.x)
                        )
               )
             )
           ) %>%
      dplyr::mutate(
                 entropy = purrr::map2(
                                      lqt, lpst,
                                      ~ calculateEntropyTwoCt(.x, .y, lambda)
                                  )
             ) %>%
      dplyr::select(var, ct, entropy, median.change.num) %>% unnest()
}


##' Calculate entropy of two change points
##'
##' This is used for determin significant change variates
##' @title calculateEntropyTwoCt
##' @param lqt Numeric vector, Each value represents Yt:n probability
##' @param lpst Numeric matrix, Log probability for each intervals
##' @param lambda Numeric vector, Each value represents Yt:n probability
##' @return entropy Numeric, entropy of two change points
##' @author Yasuhiro Kojima
calculateEntropyTwoCt <- function(lqt, lpst, lambda){
    entropy <- 0
    for(ct1 in seq(length(lqt) - 1)){
        for(ct2 in seq(ct1 + 1, length(lqt))){
            if(ct2 < length(lqt)){
                lprob <- lpst[1, ct1] + lpst[ct1 + 1, ct2] + lqt[ct2 + 1] +
                    2 * log(lambda) + (ct2 - 2) * log(1 - lambda)
            }else{
                lprob <- lpst[1, ct1] + lpst[ct1 + 1, ct2] +
                    log(lambda) + (ct2 - 1) * log(1 - lambda)
            }
            entropy <- exp(lprob) * lprob
        }
    }
    return(entropy)
}
