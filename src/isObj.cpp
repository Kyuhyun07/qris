#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export(rng = FALSE)]]
arma::mat isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H,
		arma::vec I, arma::vec logT, double Q) {
  arma::mat m1 = X % repmat(I, 1, X.n_cols);
  arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
  return m1.t() * m2;
}
