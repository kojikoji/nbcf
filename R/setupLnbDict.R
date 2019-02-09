
##' Set up \code{\link{LnbDict}} instance
##'
##' This function is responsible for generating instance of \code{\link{LnbDict}}. This involves calculate negative binomial probability for grid points of count and p. Futhermore, this calculate index conversion for count and p (see details in \code{\link{LnbDict}}. This also calculate log prior values for \code{p.vec}.
##' @title setupLnbDict
##' @param count.max dgCMatrix, Max values of \code{count.mat} in \code{link{Nbcd}}
##' @param mean.bias.vec Numeric vecotor, Relative mean bias for each observation 
##' @param r Numeric, Size parameters of negative binomial distribution. 
##' @param alpha Numeric, Parameter fo beta distribution, which is prior for p of negative binomial distribution
##' @param beta Numeric, Parameter fo beta distribution, which is prior for p of negative binomial distribution
##' @param count.res Integer, Approximation resolution of count
##' @param p.res Integer, Approximation resolution of parameter p 
##' @return lnb.dict Instance of class \code{\link{LnbDict}}
##' @seealso [LnbDict]
##' @author Yasuhiro Kojima

setupLnbDict <- function(count.max, mean.bias.vec, r, alpha, beta, count.res, p.res){
  count.vec <- calculateCountVec(count.max, count.res)
  count.dict <- calculateCountDict(count.vec)
  p.vec <- calculatePvec(p.res)
  p.width.vec <- calculatePwidthVec(p.vec)
  p.dict <- calculatePdict(p.vec, mean.bias.vec)
  lnb.values <- calculateLnbValues(count.vec, p.vec, r)
  lprior.values <- log(dbeta(p.vec, alpha, beta))
  lnb.dict <- new("LnbDict", values=lnb.values, lprior.values=lprior.values, count.dict=count.dict, p.dict=p.dict, p.width.vec=p.width.vec)
  return(lnb.dict)
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
    ## this return 0 ~ continuous.upper-2, max.val
    count.vec <- c(seq(0, count.res -2), count.max)
    warning(paste("Count values over", count.res, "will be approximated to", count.res,"or", count.max))
  }   
  return(count.vec)
}

##' Calculate count dictionary from row value to approximated points
##'
##' count dictionary is conversion vector from raw count to approximated count
##' @title calculateCountDict
##' @param count.vec Interger vector, Selected count grid points
##' @return count.dict
##' @author Yasuhiro Kojima
calculateCountDict <- function(count.vec){
  count.dict <- rcpp_calculate_count_dict(as.integer(count.vec))
  return(count.dict)
}

##' Calculation for grid points of p
##'
##' These grids points are concentrated on 0 or 1. Use quadratic grids (Ferrer et al. 2016)
##' @title calculatePvec
##' @param p.res Integer, Approximation resolution of p
##' @return p.vec Numeric vector, Grid points for p
##' @author Yasuhiro Kojima
calculatePvec <- function(p.res){
  p.vec <- vector()
  p.vec[1] <- 0
  for(i in seq(p.res-1)){
    x = i/p.res
    p.vec[i+1] <- p.vec[i] + x*(1-x)
  }
  p.vec <- p.vec/p.vec[length(p.vec)]
  return(p.vec)
}

##' Calculation for grid width of p
##'
##' Calculation is based on \code{p.vec}
##' @title calculatePwidthVec
##' @param p.vec Numeric vector, Grid points for p
##' @return p.width.vec grid width of p
##' @author Yasuhiro Kojima
calculatePwidthVec <- function(p.vec){
  grid.num <- length(p.vec)
  p.width.vec <- vector()
  p.width.vec[1] <- (p.vec[2] - p.vec[1])/2
  p.width.vec[grid.num] <- (p.vec[grid.num] - p.vec[grid.num - 1])/2
  for(i in seq(2, grid.num - 1)){
    p.width.vec[i] <- (p.vec[i+1] - p.vec[i])/2 + (p.vec[i] - p.vec[i-1])/2
  }
  return(p.width.vec)
}

##' Calculate p dictionary from each observation grid to grids corrected for mean.bias.vec
##'
##' We have to calculate negative binomial probability to p corrected for mean.bias.vec Hence, we need to use conversion before after correction.
##' @title calculatePvec
##' @param p.res Integer, Approximation resolution of p
##' @return p.dict Numeric matrix, conversion index form each grid point in each observation to reference grid points. each row represents grid point. each column represents each observation.
##' @author Yasuhiro Kojima
calculatePdict <- function(p.vec, mean.bias.vec){
  mean.bias.vec <- mean.bias.vec/mean(mean.bias.vec)
  p.dict <- rcpp_calculate_p_dict(p.vec, mean.bias.vec)
  return(p.dict)
}

##' Calculate the probability value of log negative binomial distribution
##'
##' The values are calculated for each count (row) and p (column)
##' @title calculateLnbValues
##' @param count.vec Interger vector, Selected count grid points
##' @param p.vec Numeric vector, Selected p grid points
##' @param r Numeric, Size parameters of negative binomial distribution. 
##' @return lnb.mat Numeric matrix, log negative binomial probability. Columns and rows correspond to grid points of count and p
##' @author Yasuhiro Kojima
calculateLnbValues <- function(count.vec, p.vec, r){
  lnb.mat <- rcpp_calculate_lnb_values(count.vec, p.vec, r)
  return(lnb.mat)
}
