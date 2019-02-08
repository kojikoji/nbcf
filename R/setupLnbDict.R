
##' Set up \code{\link{LnbDict}} instance
##'
##' This function is responsible for generating instance of \code{\link{LnbDict}}. This involves calculate negative binomial probability for grid points of count and p. Futhermore, this calculate index conversion for count and p (see details in \code{\link{LnbDict}}.
##' @title setupLnbDict
##' @param count.max dgCMatrix, Max values of \code{count.mat} in \code{link{Nbcd}}
##' @param mean.bias.vec Numeric vecotor, Relative mean bias for each observation 
##' @param p.res Integer, Approximation resolution of parameter p 
##' @param count.res Integer, Approximation resolution of count
##' @return lnbdict Instance of class \code{\link{LnbDict}}
##' @seealso [LnbDict]
##' @author Yasuhiro Kojima

setupLnbDict <- function(count.max, p.res, count.res, mean.bias.vec){
  count.vec <- calculateCountVec(count.max, count.res)
  count.dict <- calculateCountDict(count.vec)
  p.vec <- calculatePvec(p.res)
  p.dict <- calculatePdict(p.vec, mean.bias)
  lnb.values <- calculateLnbValue(p.vec, count.vec)
  lnbdict <- new("LnbDict", values=lnb.values, count.dict=count.dict, p.dict=p.dict)
  return(lnbdict)
}

##' Chose count values for approximation
##'
##' Chose \code{count.res} count values from 0 ~ \code{max(count.mat)} for approximation. This chose all values up to min(continuos.upper, max(count.mat)) and chose other values based on bining of log(count).
##' @title calculateCountVec
##' @param count.max Integer, Max values of \code{count.mat} in \code{link{Nbcd}}
##' @param count.res Integer, Approximation resolution of count
##' @return count.vec Interger vector, Selected count grid points
##' @author Yasuhiro Kojima

calculateCountVec <- function(count.max, count.res, continuous.upper=300){
  if(count.max < count.res) return(seq(0, count.max))
  # continuous value are selected until upper
  continuous.upper <- min(count.max, continuous.upper)
  continuous.vec <- seq(0, continuous.upper)
  if(continuous.upper < count.res){
    ## over the threshold we select values from continuous + 1 to max.count
    ## at regular interval based on log scale
    large.grid.num <- count.res - length(continuous.vec)
    log.max.val <- log(count.max)
    log.min.val <- log(continuous.upper + 1)
    log.large.vec <- seq(log.min.val, log.max.val, length.out = large.grid.num)
    large.vec <- round(exp(log.large.vec))
    ## return values are concatenated above continuous and selected values
    count.vec <- c(continuous.vec, large.vec)
  }else{
    ## if continuous.upper is larger than count.res,
    ## this return 0 ~ continuous.upper
    count.vec <- seq(0, count.res -1)
    warning(paste("Count values over", count.res, "will be approximated to", count.res))
  }   
  return(count.vec)
}
