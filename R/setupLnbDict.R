
##' Set up \code{\link{LnbDict}} instance
##'
##' This function is responsible for generating instance of \code{\link{LnbDict}}. This involves calculate negative binomial probability for grid points of count and p. Futhermore, this calculate index conversion for count and p (see details in \code{\link{LnbDict}}.
##' @title setupLnbDict
##' @param count.max Max values of \code{count.mat} in \code{link{Nbcd}}
##' @param mean.bias.vec Relative mean bias for each observation 
##' @param p.res Approximation resolution of parameter p 
##' @param count.res Approximation resolution of count
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
