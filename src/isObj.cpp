#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export(rng = FALSE)]]
arma::mat isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H,
		arma::vec I, arma::vec logT, double Q) {
  // arma::mat m1 = X % repmat(I, 1, X.n_cols);
  // arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
  arma::mat m1 = X;
  m1.each_col() %= I;
  arma::mat m2 = normcdf((X * b - logT) / sqrt(sum(X % (X * H), 1))) % W - Q;
  return m1.t() * m2;
}

//' @noRd
// [[Rcpp::export(rng = FALSE)]]
Rcpp::List isObjL(arma::vec b, arma::mat X, arma::vec W, arma::mat H,
  		 arma::vec I, arma::vec logT, double Q) {
  // arma::mat m1 = X % repmat(I, 1, X.n_cols);
  // arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
  arma::mat m1 = X;
  m1.each_col() %= I;
  Rcpp::List out(2);
  out(0) = m1;
  out(1) = normcdf((X * b - logT) / sqrt(sum(X % (X * H), 1)));
  out.names() = Rcpp::CharacterVector::create("m1", "m2");
  return(out); 
}
