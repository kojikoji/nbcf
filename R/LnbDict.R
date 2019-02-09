## LnbDict class
## Yasuhiro Kojima

##' Data stroge for log negative binomial probability
##'
##' This class is a storage for precomputated log negative binomial probability. This also has the index conversion matrix and vector for parameter p and count.
##'
##' @title Lnbdict
##' @slot values Numeric Matrix, This storage log negative binomial probability. Each column and row represent each grid of count and p.
##' @slot lprior.values Numeric vector, log prior probability values for each gridpoints
##' @slot count.dict This is a sparse vector each elements have the row index of \code{value} cressponding the count represented by index of this vector
##' @slot p.dict This is a dense matrix. Each element has a column index of \code{value} after correction against \code{\link[Nbcf]{mean.bias.vec}}. Each row represents index of p before correction. Each column represents observation index.
##' @author Yasuhiro Kojima
##' 
##' @import Matrix

setClass(
  "LnbDict",
  representation(
    values="matrix",
    lprior.values="vector",
    count.dict="vector",
    p.dict="matrix")
)
