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
  if(length(ct.vec) > 0){
    ## add ct vec to begin and end of time
    extended.ct.vec <- c(0, ct.vec, ncol(lpst.list[[1]]))
    change.var.df <- tibble()
    for(var in used.vars){
      for(ct.idx in seq(length(ct.vec))){
        ct  <-  ct.vec[ct.idx]
        previous.ct  <-  extended.ct.vec[ct.idx - 1 + 1]
        following.ct  <-  extended.ct.vec[ct.idx + 1 + 1]
        lpst  <-  lpst.list[[var]][(previous.ct + 1):following.ct,
          (previous.ct + 1):following.ct]
        lqt  <-  calculateLqt(lpst, lambda)
        sim.ct.list  <-  perfectSimulation(lqt, lpst, lambda)
        median.change.num  <-  median(
          unlist(purrr::map(
                          sim.ct.list,
                          ~ length(.x)
                        )
                 )
        )
        entropy <- calculateEntropyTwoCt(lqt, lpst, lambda)
        change.var.df <- rbind(change.var.df,
                               tibble(ct = ct,
                                      var = var,
                                      median.change.num = median.change.num,
                                      entropy = entropy
                                      )
                               )
      }
    }
  }else{
    change.var.df <- tibble()
  }
  return(change.var.df)
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
    regulalize <- function(vec) vec/sum(vec)
    regExp <- function(vec) regulalize(exp(vec - max(vec)))
    lprob.vec <- vector()
    for(ct1 in seq(length(lqt) - 1)){
        for(ct2 in seq(ct1 + 1, length(lqt))){
            if(ct2 < length(lqt)){
                lprob <- lpst[1, ct1] + lpst[ct1 + 1, ct2] + lqt[ct2 + 1] +
                    2 * log(lambda) + (ct2 - 2) * log(1 - lambda)
            }else{
                lprob <- lpst[1, ct1] + lpst[ct1 + 1, ct2] +
                    log(lambda) + (ct2 - 1) * log(1 - lambda)
            }
            lprob.vec <- c(lprob.vec, lprob)
        }
    }
    prob.vec <- regExp(lprob.vec)
    entropy <- sum(- prob.vec * log(prob.vec + 1.0e-200))
    return(entropy)
}
