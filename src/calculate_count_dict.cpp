#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

// Calculate count dictionary
//
// count dictionary is conversion vector from raw count to approximated count
// @title calculateCountDict
// @param count.vec Interger vector, Selected count grid points
// @return count.dict
// @author Yasuhiro Kojima

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


Eigen::VectorXi rcpp_calculate_count_dict(const Eigen::Map<Eigen::VectorXi> count_vec){
  int count_max = count_vec.maxCoeff();
  Eigen::VectorXi count_dict(count_max+1);
  // initialize grid bound
  // We select two grids between which each original values are
  int grid_lower = count_vec(0);
  int grid_upper = count_vec(1);
  int grid_upper_index = 1;
  for(int i = 0; i < count_max+1; i++){
    if(i > grid_upper){
      // Shift grid
      grid_lower = grid_upper;
      grid_upper = count_vec(++grid_upper_index);
    }
    // Chose lower or upper grid
    if(i - grid_lower < grid_upper - i){
      count_dict(i) = grid_lower;
    }else{
      count_dict(i) = grid_upper;
    }
  }
  return(count_dict);
}
