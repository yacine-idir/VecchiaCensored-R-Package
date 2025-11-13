// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export("precision_decomp_sparse")]]
arma::sp_mat precision_decomp_sparse(List cov_par,
                                     const arma::mat& X,
                                     const IntegerMatrix& NN,
                                     const arma::mat& locs,
                                     int start = 1) {  // Added start parameter with default = 1
  const int n = locs.n_rows;
  const int p = X.n_cols;
  
  // Convert R's 1-based indexing to C++'s 0-based indexing
  const int start_idx = start - 1;
  
  // parse covariance parameters
  std::vector<double> sigma2(p), phi(p);
  for (int k = 0; k < p; ++k) {
    NumericVector vk = cov_par[k];
    sigma2[k] = vk["sigma_sq"];
    phi[k]    = vk["phi"];
  }
  double nugget = as<double>(cov_par["nugget"]);
  
  arma::sp_mat U(n, n);
  
  // handle first point (only if start_idx == 0)
  if (start_idx == 0) {
    double c0 = 0.0;
    for (int k = 0; k < p; ++k) {
      c0 += X(0, k) * X(0, k) * sigma2[k];
    }
    c0 += nugget;
    U(0, 0) = 1.0 / std::sqrt(c0);
  }
  // If start_idx > 0, row 0 remains all zeros
  
  // temporary storage
  arma::mat Cmat;
  arma::mat cov_nn;
  arma::rowvec cov_i_nn;
  arma::mat L;
  
  // parallel over rows, but only starting from max(1, start_idx)
  int loop_start = std::max(1, start_idx);
  
#pragma omp parallel for private(Cmat, cov_nn, cov_i_nn, L)
  for (int i = loop_start; i < n; ++i) {  // Start from loop_start instead of 1
    std::vector<int> neigh;
    for (int jj = 0; jj < NN.ncol(); ++jj) {
      int r = NN(i, jj);
      if (r == NA_INTEGER) break;
      neigh.push_back(r - 1);
    }
    int M = neigh.size();
    if (M < 2) {
      U(i, i) = 1.0 / std::sqrt(nugget);
      continue;
    }
    
    // build neighborhood covariance matrix
    Cmat.set_size(M, M);
    for (int u = 0; u < M; ++u) {
      const arma::rowvec loc_u = locs.row(neigh[u]);
      for (int v = u; v < M; ++v) {
        const arma::rowvec loc_v = locs.row(neigh[v]);
        double cc = 0.0;
        double dist = std::sqrt(arma::dot(loc_u - loc_v, loc_u - loc_v));
        for (int k = 0; k < p; ++k) {
          cc += X(neigh[u], k) * X(neigh[v], k) * std::exp(-dist / phi[k]) * sigma2[k];
        }
        if (u == v) cc += nugget;
        Cmat(u, v) = Cmat(v, u) = cc;
      }
    }
    
    cov_nn = Cmat.submat(1, 1, M-1, M-1);
    cov_i_nn = Cmat.submat(0, 1, 0, M-1);
    
    arma::chol(L, cov_nn, "lower");
    arma::vec y = arma::solve(arma::trimatl(L), cov_i_nn.t());
    arma::vec x = arma::solve(arma::trimatu(L.t()), y);
    double d0 = Cmat(0, 0) - arma::dot(cov_i_nn, x);
    double sd = 1.0 / std::sqrt(d0);
    
    U(i, i) = sd;
    for (int k = 0; k < M-1; ++k) {
      U(i, neigh[k+1]) = -x[k] * sd;
    }
  }
  
  U.sync();
  return U;
}