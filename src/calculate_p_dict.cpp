#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
// Calculate p_vec after correction of mean_bias
VectorXd fix_p_vec(VectorXd p_vec, double mean_bias){
  VectorXd fixed_p_vec = p_vec.unaryExpr([=](double p){return(mean_bias*p/(1+(mean_bias-1)*p));});
  return(fixed_p_vec);
}

// Get approximated index of each of original_vec in approximate_vec.
// Note: This function assume both vectors sorted ascending
// [[Rcpp::export]]

Eigen::VectorXi get_approximate_index(Eigen::VectorXd original_vec, Eigen::VectorXd approximate_vec){
  Eigen::VectorXi approximate_index(original_vec.size());
  // initialize grid bound
  // We select two grids between which each original values are
  double grid_lower = approximate_vec(0);
  double grid_upper = approximate_vec(1);
  int grid_upper_index = 1;
  for(int i = 0; i < original_vec.size(); i++){
    // Once a value of original exceeds max of approximated,
    // left values (bigger ones) are all approximated max of approximated
    if(original_vec(i) > approximate_vec.maxCoeff()){
      approximate_index(i) = approximate_vec.size() - 1;
      continue;
    }
    if(original_vec(i) > grid_upper){
      // Shift grid
      while(original_vec(i) > grid_upper && grid_upper_index < approximate_vec.size() - 1){
	grid_lower = grid_upper;
	grid_upper = approximate_vec(++grid_upper_index);
      }
    }
    // Chose lower or upper grid
    if(original_vec(i) - grid_lower < grid_upper - original_vec(i)){
      approximate_index(i) = grid_upper_index -1;
    }else{
      approximate_index(i) = grid_upper_index;
    }
  }
  return(approximate_index);
}


// Calculate p dictionary from each observation grid to grids corrected for mean.bias.vec
//
// We have to calculate negative binomial probability to p corrected for mean.bias.vec Hence, we need to use conversion before after correction.
// @title calculatePvec
// @param p.res Integer, Approximation resolution of p
// @return p.dict Numeric matrix, conversion index form each grid point in each observation to reference grid points. each row represents grid point. each column represents each observation.
// @author Yasuhiro Kojima

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXi rcpp_calculate_p_dict(const Eigen::Map<Eigen::VectorXd> p_vec, const Eigen::Map<Eigen::VectorXd> mean_bias_vec){
  int obs_num = mean_bias_vec.size();
  int grid_num = p_vec.size();
  Eigen::MatrixXi p_dict(grid_num, obs_num);
  for(int i = 0; i < obs_num; i++){
    Eigen::VectorXd fixed_p_vec = fix_p_vec(p_vec, mean_bias_vec(i));
    p_dict.col(i) = get_approximate_index(fixed_p_vec, p_vec);
  }
  return(p_dict);
}


