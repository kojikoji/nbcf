#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>

// Calculate the probability value of log negative binomial distribution
//
// The values are calculated for each count (column) and p (row)
// @title calculateLnbValues
// @param count.vec Interger vector, Selected count grid points
// @param p.vec Interger vector, Selected p grid points
// @return lnb.mat Numeric matrix, log negative binomial probability. Columns and rows correspond to grid points of count and p 
// @author Yasuhiro Kojima

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd rcpp_calculate_lnb_values(Eigen::VectorXi count_vec, Eigen::VectorXd p_vec, double r){
  Eigen::MatrixXd lnb_values(p_vec.size(), count_vec.size());
  for(int i = 0; i < count_vec.size(); i++){
    for(int j = 0; j < p_vec.size(); j++){
      double k = count_vec(i);
      double p = p_vec(j);
      // When p ==  0 or 1, avoid emergence of infinity
      double logp = p == 0 ? -1.0e30 : log(p);
      double log1mp = (1 - p) == 0 ? -1.0e30 : log(1 - p);
      lnb_values(j, i) = lgamma(r + k) - lgamma(r) - lgamma(k+1) + k * logp + r * log1mp;
    }
  }
  return(lnb_values);
}
