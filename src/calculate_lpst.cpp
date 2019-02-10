#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;

// p_grid * t_grid lnb value matrix. lnb values within same t_grid is summarized
// [[Rcpp::export]]
Eigen::MatrixXd calculate_pt_grid_lp(const Eigen::VectorXd &count_sp_vec,
			      const Eigen::MatrixXd &grid_lnb_mat,
			      const Eigen::VectorXi &count_dict,
			      const Eigen::MatrixXi &p_dict,
			      const Eigen::VectorXi &t_grids){
  Eigen::MatrixXd pt_grid_lp = MatrixXd::Zero(p_dict.rows(), t_grids.size() - 1);
  int upper_t_grid_id = 1;
  for(int oid = 0; oid < count_sp_vec.size(); oid++){
    // when oid is out of current grid interval, change grid interval
    if(oid > t_grids(upper_t_grid_id)){
      upper_t_grid_id++;
    }
    for(int p_id = 0; p_id < p_dict.rows(); p_id++){
      // add the lnb value to current grid
      int approx_count_id = count_dict(count_sp_vec(oid));
      pt_grid_lp(p_id, upper_t_grid_id - 1) += grid_lnb_mat(p_id, approx_count_id);
    }
  }
  return(pt_grid_lp);
}

// Accumulate the values of mat column wise
// [[Rcpp::export]]
Eigen::MatrixXd accumulate_colwise(const Eigen::MatrixXd &mat){
  Eigen::MatrixXd accumulated_mat = Eigen::MatrixXd(mat.rows(), mat.cols());
  accumulated_mat.col(0) = mat.col(0);
  for(int coli = 1; coli < mat.cols(); coli++){
    accumulated_mat.col(coli) = accumulated_mat.col(coli - 1) + mat.col(coli);
  }
  return(accumulated_mat);
}


// summation is based on w_vec
// [[Rcpp::export]]
double log_sum_exp(const Eigen::VectorXd &vec, const Eigen::VectorXd &w_vec){
  double max_val = vec.maxCoeff();
  double val = log(((vec.array() - max_val).exp() * w_vec.array()).sum()) + max_val;
  return(val);
}

// summation is based on w_vec
// [[Rcpp::export]]
Eigen::VectorXd colwise_log_sum_exp(const Eigen::MatrixXd &mat, const Eigen::VectorXd &w_vec){
  Eigen::VectorXd applied_vec(mat.cols());
  for(int coli = 0; coli < mat.cols(); coli++){
    applied_vec(coli) = log_sum_exp(mat.col(coli), w_vec);
  }
  return(applied_vec);
}


// calculate lpst for each variate
// [[Rcpp::export]]
Eigen::MatrixXd calculate_lpst_g(const Eigen::VectorXd &count_sp_vec,
				 const Eigen::MatrixXd &grid_lnb_mat,
				 const Eigen::VectorXd &lprior_vec,
				 const Eigen::VectorXi &count_dict,
				 const Eigen::MatrixXi &p_dict,
				 const Eigen::VectorXd &p_width_vec,
				 const Eigen::VectorXi &t_grids){
  Eigen::MatrixXd pt_grid_lp = calculate_pt_grid_lp(count_sp_vec, grid_lnb_mat, count_dict, p_dict, t_grids);
  // Accumulation make pt_grid_lp summed log probability from 0 to each column index
  Eigen::MatrixXd p_grid_lpsn = accumulate_colwise(pt_grid_lp);
  Eigen::MatrixXd lp_mat = Eigen::MatrixXd::Zero(t_grids.size() - 1,  t_grids.size() - 1);
  for(int s = 0; s < t_grids.size() - 1; s++){
    Eigen::MatrixXd pr_p_grid_lpsn = p_grid_lpsn.colwise() + lprior_vec; 
    lp_mat.row(s) = colwise_log_sum_exp(pr_p_grid_lpsn, p_width_vec);
    p_grid_lpsn = p_grid_lpsn.colwise() - pt_grid_lp.col(s);
  }
  return(lp_mat);
}

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_calculate_lpst(const Eigen::SparseMatrix<double> &count_sp_mat,
				    const Eigen::MatrixXd &grid_lnb_mat,
				    const Eigen::VectorXd &lprior_vec,
				    const Eigen::VectorXi &count_dict,
				    const Eigen::MatrixXi &p_dict,
				    const Eigen::VectorXd &p_width_vec,
				    const Eigen::VectorXi &t_grids){
  Eigen::MatrixXd lpst = Eigen::MatrixXd::Zero(t_grids.size() - 1,  t_grids.size() - 1);
  for(int g = 0; g < count_sp_mat.rows(); g++){
    lpst += calculate_lpst_g(count_sp_mat.row(g), grid_lnb_mat, lprior_vec, count_dict, p_dict, p_width_vec, t_grids);
  }
  return(lpst);
}
