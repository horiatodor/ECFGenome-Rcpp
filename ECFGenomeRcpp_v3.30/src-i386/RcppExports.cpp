// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// do_operon_Rcpp
NumericMatrix do_operon_Rcpp(NumericVector all_scores, List operon_two, IntegerVector matches_to_all_scores);
RcppExport SEXP _ECFGenomeRcpp_do_operon_Rcpp(SEXP all_scoresSEXP, SEXP operon_twoSEXP, SEXP matches_to_all_scoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type all_scores(all_scoresSEXP);
    Rcpp::traits::input_parameter< List >::type operon_two(operon_twoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type matches_to_all_scores(matches_to_all_scoresSEXP);
    rcpp_result_gen = Rcpp::wrap(do_operon_Rcpp(all_scores, operon_two, matches_to_all_scores));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _ECFGenomeRcpp_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _ECFGenomeRcpp_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _ECFGenomeRcpp_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _ECFGenomeRcpp_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _ECFGenomeRcpp_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// scoreseq_plus_single2
List scoreseq_plus_single2(NumericMatrix pwm35, ListOf<CharacterVector> sequences);
RcppExport SEXP _ECFGenomeRcpp_scoreseq_plus_single2(SEXP pwm35SEXP, SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type pwm35(pwm35SEXP);
    Rcpp::traits::input_parameter< ListOf<CharacterVector> >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreseq_plus_single2(pwm35, sequences));
    return rcpp_result_gen;
END_RCPP
}
// scoreseq_plus_combinatoric2
List scoreseq_plus_combinatoric2(NumericMatrix pwm35, NumericMatrix pwm10, ListOf<CharacterVector> sequences, IntegerVector spacer);
RcppExport SEXP _ECFGenomeRcpp_scoreseq_plus_combinatoric2(SEXP pwm35SEXP, SEXP pwm10SEXP, SEXP sequencesSEXP, SEXP spacerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type pwm35(pwm35SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pwm10(pwm10SEXP);
    Rcpp::traits::input_parameter< ListOf<CharacterVector> >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type spacer(spacerSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreseq_plus_combinatoric2(pwm35, pwm10, sequences, spacer));
    return rcpp_result_gen;
END_RCPP
}
// sum_which_greater
NumericVector sum_which_greater(NumericVector of_int, NumericVector to_sum, NumericVector avector);
RcppExport SEXP _ECFGenomeRcpp_sum_which_greater(SEXP of_intSEXP, SEXP to_sumSEXP, SEXP avectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type of_int(of_intSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type to_sum(to_sumSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type avector(avectorSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_which_greater(of_int, to_sum, avector));
    return rcpp_result_gen;
END_RCPP
}
// sum_which_greater_or_equal
NumericVector sum_which_greater_or_equal(NumericVector x, NumericVector y, NumericVector z);
RcppExport SEXP _ECFGenomeRcpp_sum_which_greater_or_equal(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_which_greater_or_equal(x, y, z));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ECFGenomeRcpp_do_operon_Rcpp", (DL_FUNC) &_ECFGenomeRcpp_do_operon_Rcpp, 3},
    {"_ECFGenomeRcpp_rcpp_hello_world", (DL_FUNC) &_ECFGenomeRcpp_rcpp_hello_world, 0},
    {"_ECFGenomeRcpp_rcpparma_hello_world", (DL_FUNC) &_ECFGenomeRcpp_rcpparma_hello_world, 0},
    {"_ECFGenomeRcpp_rcpparma_outerproduct", (DL_FUNC) &_ECFGenomeRcpp_rcpparma_outerproduct, 1},
    {"_ECFGenomeRcpp_rcpparma_innerproduct", (DL_FUNC) &_ECFGenomeRcpp_rcpparma_innerproduct, 1},
    {"_ECFGenomeRcpp_rcpparma_bothproducts", (DL_FUNC) &_ECFGenomeRcpp_rcpparma_bothproducts, 1},
    {"_ECFGenomeRcpp_scoreseq_plus_single2", (DL_FUNC) &_ECFGenomeRcpp_scoreseq_plus_single2, 2},
    {"_ECFGenomeRcpp_scoreseq_plus_combinatoric2", (DL_FUNC) &_ECFGenomeRcpp_scoreseq_plus_combinatoric2, 4},
    {"_ECFGenomeRcpp_sum_which_greater", (DL_FUNC) &_ECFGenomeRcpp_sum_which_greater, 3},
    {"_ECFGenomeRcpp_sum_which_greater_or_equal", (DL_FUNC) &_ECFGenomeRcpp_sum_which_greater_or_equal, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ECFGenomeRcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
