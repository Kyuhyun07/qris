#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export(rng = FALSE)]]
arma::mat Amat(arma::vec b, arma::mat X, arma::vec W_star, arma::mat H,
	       arma::vec I,arma::vec logT, double Q) {
  arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(W_star, 1, X.n_cols);
  arma::mat m2 =  (-X / repmat(sqrt(diagvec(X * H * X.t())), 1, X.n_cols)) %
    repmat(normpdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))), 1, X.n_cols);
  return m1.t() * m2;
}
