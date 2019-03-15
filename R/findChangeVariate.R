##' Find change variates at each change points
##'
##' This function calculate the signicance of  change points for each variate variates.
##' @title findChangeVariate
##' @param lpst.list list of Numeric matrix, list of log probability for each intervals for each variate
##' @param used.vars vector, names of variates used in change point finding
##' @param ct.vec Numeric vector, change point vector
##' @return change.var.df tibble, contains variate names, change points and bayes factors
##' @author Yasuhiro Kojima
##' @import tidyverse
findChangeVariate <- function(lpst.list, used.vars, ct.vec, lambda = 0.5){
  if(length(ct.vec) > 0){
    ## add ct vec to begin and end of time
    extended.ct.vec <- c(1, ct.vec, ncol(lpst.list[[1]]))
    ##:ess-bp-start::browser@nil:##
    purrr::map_dfr(
             cross(list(var = used.vars, ct.idx = seq(length(ct.vec)))),
             function(var.list){
               var <- var.list[["var"]]
               ct.idx <- var.list[["ct.idx"]]
               ct  <-  ct.vec[ct.idx]
               previous.ct  <-  extended.ct.vec[ct.idx - 1 + 1]
               following.ct  <-  extended.ct.vec[ct.idx + 1 + 1]
               joined.lpst <- lpst.list[[var]][previous.ct, following.ct]
               divided.lpst <- lpst.list[[var]][previous.ct, ct] + lpst.list[[var]][(ct+1), following.ct]
               bayes.factor <- divided.lpst - joined.lpst
               this.length <- following.ct - previous.ct
               total.length <- ncol(lpst.list[[var]])
               var.ratio <- (-bayes.factor/this.length)/(lpst.list[[var]][1,total.length]/total.length)
               ##:ess-bp-start::browser@nil:##
tibble(var = var, change.point = ct, bayes.factor = bayes.factor, var.ratio = var.ratio)
             }
           )
  }else{
    tibble()
  }
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
