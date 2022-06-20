// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Amat
arma::mat Amat(arma::vec b, arma::mat X, arma::vec W_star, arma::mat H, arma::vec I, arma::vec logT, double Q);
RcppExport SEXP _qris_Amat(SEXP bSEXP, SEXP XSEXP, SEXP W_starSEXP, SEXP HSEXP, SEXP ISEXP, SEXP logTSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W_star(W_starSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logT(logTSEXP);
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(Amat(b, X, W_star, H, I, logT, Q));
    return rcpp_result_gen;
END_RCPP
}
// ghat
arma::vec ghat(arma::vec Time, arma::vec censor, arma::vec wgt);
RcppExport SEXP _qris_ghat(SEXP TimeSEXP, SEXP censorSEXP, SEXP wgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type censor(censorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wgt(wgtSEXP);
    rcpp_result_gen = Rcpp::wrap(ghat(Time, censor, wgt));
    return rcpp_result_gen;
END_RCPP
}
// isObj
arma::mat isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec I, arma::vec logT, double Q);
RcppExport SEXP _qris_isObj(SEXP bSEXP, SEXP XSEXP, SEXP WSEXP, SEXP HSEXP, SEXP ISEXP, SEXP logTSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logT(logTSEXP);
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(isObj(b, X, W, H, I, logT, Q));
    return rcpp_result_gen;
END_RCPP
}
// isObjL
Rcpp::List isObjL(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec I, arma::vec logT, double Q);
RcppExport SEXP _qris_isObjL(SEXP bSEXP, SEXP XSEXP, SEXP WSEXP, SEXP HSEXP, SEXP ISEXP, SEXP logTSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logT(logTSEXP);
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(isObjL(b, X, W, H, I, logT, Q));
    return rcpp_result_gen;
END_RCPP
}
// rev_isObj
arma::mat rev_isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec E, arma::vec I, arma::vec logT, double Q);
RcppExport SEXP _qris_rev_isObj(SEXP bSEXP, SEXP XSEXP, SEXP WSEXP, SEXP HSEXP, SEXP ESEXP, SEXP ISEXP, SEXP logTSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logT(logTSEXP);
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(rev_isObj(b, X, W, H, E, I, logT, Q));
    return rcpp_result_gen;
END_RCPP
}
