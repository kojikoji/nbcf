##' Conduct perfectsimulation for change point
##'
##' This function samples change points using perfectsimulation, which is a machine learning techunique.
##' @title perfectSimulation
##' @param lqt Numeric vector, Each value represents Yt:n probability 
##' @param lpst Numeric matrix, Log probability for each intervals
##' @param lambda Numeric vector, Each value represents Yt:n probability
##' @param sim.iter Integer, Iteration of simulation
##' @return change.points.list List, list of numeric vector of change points at each simulation
##' @author Yasuhiro Kojima
perfectSimulation <- function(lqt, lpst, lambda, sim.iter=100){
  transition.prob.list <- calculateCtTransitionProb(lqt, lpst, lambda)
  change.points.list <- list()
  t.max <- length(lqt)
  for(i in seq(sim.iter)){
    ct.vec <- vector()
    current.ct <- 0
    while(current.ct < t.max){
      ## when current.ct == t.max, sample do not work
      current.ct <- ifelse(
        current.ct < t.max - 1,
        sample(seq(current.ct + 1, t.max), 1, prob=transition.prob.list[[as.character(current.ct)]]), ## as.character used because use 0 index
        t.max)
      ct.vec <- c(ct.vec, current.ct)
    }
    change.points.list[[i]] <- ct.vec
  }
  return(change.points.list)
}

##' Calculate transition probability for change points
##'
##' This calculate the probability of next change points when we assume current change points are on each time points.
##' @title calculateCtTransitionProb
##' @param lqt Numeric vector, Each value represents Yt:n probability 
##' @param lpst Numeric matrix, Log probability for each intervals
##' @param lambda Numeric vector, Each value represents Yt:n probability
##' @return transition.prob.list List, each element represents vector of probability of next change point. index represents previous change points
##' @author Yasuhiro Kojima
calculateCtTransitionProb <- function(lqt, lpst, lambda){
  pre.ct <- 0
  t.max <- length(lqt)
  transition.prob.list <- list()
  for(pre.ct in seq(0, t.max-1)){
    ## when pre.t = t.max -1, there is no futhre change
    if(pre.ct < t.max -1){
      lprob.vec <- unlist(purrr::map((pre.ct+1):(t.max-1),
                                     ~ lpst[(pre.ct+1),.x] + (.x - pre.ct) * log(1-lambda) +
                                       log(lambda) + lqt[.x + 1]))
    }else{
      lprob.vec <- c()
    }
    non.change.lprob <- lpst[(pre.ct + 1), t.max] + (t.max - pre.ct) * log(1 - lambda)
    lprob.vec <- c(lprob.vec, non.change.lprob)
    regulalize <- function(vec) vec/sum(vec)
    regExp <- function(vec) regulalize(exp(vec - max(vec)))
    prob.vec <- regExp(lprob.vec)
    ## as.character used because use 0 index
    transition.prob.list[[as.character(pre.ct)]] <- prob.vec
  }
  return(transition.prob.list)
}
