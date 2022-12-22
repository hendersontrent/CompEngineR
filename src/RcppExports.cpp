// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ac_9
double ac_9(NumericVector acfv);
RcppExport SEXP _CompEngineR_ac_9(SEXP acfvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type acfv(acfvSEXP);
    rcpp_result_gen = Rcpp::wrap(ac_9(acfv));
    return rcpp_result_gen;
END_RCPP
}
// binarize_mean
NumericVector binarize_mean(NumericVector x);
RcppExport SEXP _CompEngineR_binarize_mean(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(binarize_mean(x));
    return rcpp_result_gen;
END_RCPP
}
// embed2_incircle
double embed2_incircle(NumericVector x, NumericVector acfv, double boundary);
RcppExport SEXP _CompEngineR_embed2_incircle(SEXP xSEXP, SEXP acfvSEXP, SEXP boundarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acfv(acfvSEXP);
    Rcpp::traits::input_parameter< double >::type boundary(boundarySEXP);
    rcpp_result_gen = Rcpp::wrap(embed2_incircle(x, acfv, boundary));
    return rcpp_result_gen;
END_RCPP
}
// f_entropy
double f_entropy(NumericVector x);
RcppExport SEXP _CompEngineR_f_entropy(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f_entropy(x));
    return rcpp_result_gen;
END_RCPP
}
// firstmin_ac
int firstmin_ac(NumericVector x, NumericVector acfv);
RcppExport SEXP _CompEngineR_firstmin_ac(SEXP xSEXP, SEXP acfvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acfv(acfvSEXP);
    rcpp_result_gen = Rcpp::wrap(firstmin_ac(x, acfv));
    return rcpp_result_gen;
END_RCPP
}
// firstzero_ac
int firstzero_ac(NumericVector x, NumericVector acfv);
RcppExport SEXP _CompEngineR_firstzero_ac(SEXP xSEXP, SEXP acfvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acfv(acfvSEXP);
    rcpp_result_gen = Rcpp::wrap(firstzero_ac(x, acfv));
    return rcpp_result_gen;
END_RCPP
}
// fluctanal_prop_r1
double fluctanal_prop_r1(NumericVector x);
RcppExport SEXP _CompEngineR_fluctanal_prop_r1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fluctanal_prop_r1(x));
    return rcpp_result_gen;
END_RCPP
}
// histogram_mode
double histogram_mode(NumericVector x, int numBins);
RcppExport SEXP _CompEngineR_histogram_mode(SEXP xSEXP, SEXP numBinsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type numBins(numBinsSEXP);
    rcpp_result_gen = Rcpp::wrap(histogram_mode(x, numBins));
    return rcpp_result_gen;
END_RCPP
}
// localsimple_taures
int localsimple_taures(NumericVector x, String forecastMeth);
RcppExport SEXP _CompEngineR_localsimple_taures(SEXP xSEXP, SEXP forecastMethSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type forecastMeth(forecastMethSEXP);
    rcpp_result_gen = Rcpp::wrap(localsimple_taures(x, forecastMeth));
    return rcpp_result_gen;
END_RCPP
}
// motiftwo_entro3
double motiftwo_entro3(NumericVector x);
RcppExport SEXP _CompEngineR_motiftwo_entro3(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(motiftwo_entro3(x));
    return rcpp_result_gen;
END_RCPP
}
// outlierinclude_mdrmd
double outlierinclude_mdrmd(NumericVector x, bool zscored);
RcppExport SEXP _CompEngineR_outlierinclude_mdrmd(SEXP xSEXP, SEXP zscoredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type zscored(zscoredSEXP);
    rcpp_result_gen = Rcpp::wrap(outlierinclude_mdrmd(x, zscored));
    return rcpp_result_gen;
END_RCPP
}
// sampen_first
double sampen_first(NumericVector x);
RcppExport SEXP _CompEngineR_sampen_first(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sampen_first(x));
    return rcpp_result_gen;
END_RCPP
}
// sampenc
double sampenc(NumericVector x, int M, double r);
RcppExport SEXP _CompEngineR_sampenc(SEXP xSEXP, SEXP MSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(sampenc(x, M, r));
    return rcpp_result_gen;
END_RCPP
}
// spreadrandomlocal_meantaul
double spreadrandomlocal_meantaul(NumericVector x, int l);
RcppExport SEXP _CompEngineR_spreadrandomlocal_meantaul(SEXP xSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(spreadrandomlocal_meantaul(x, l));
    return rcpp_result_gen;
END_RCPP
}
// std1st_der
double std1st_der(NumericVector x);
RcppExport SEXP _CompEngineR_std1st_der(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(std1st_der(x));
    return rcpp_result_gen;
END_RCPP
}
// trev_num
double trev_num(NumericVector x);
RcppExport SEXP _CompEngineR_trev_num(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(trev_num(x));
    return rcpp_result_gen;
END_RCPP
}
// walker_propcross
double walker_propcross(NumericVector y);
RcppExport SEXP _CompEngineR_walker_propcross(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(walker_propcross(y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CompEngineR_ac_9", (DL_FUNC) &_CompEngineR_ac_9, 1},
    {"_CompEngineR_binarize_mean", (DL_FUNC) &_CompEngineR_binarize_mean, 1},
    {"_CompEngineR_embed2_incircle", (DL_FUNC) &_CompEngineR_embed2_incircle, 3},
    {"_CompEngineR_f_entropy", (DL_FUNC) &_CompEngineR_f_entropy, 1},
    {"_CompEngineR_firstmin_ac", (DL_FUNC) &_CompEngineR_firstmin_ac, 2},
    {"_CompEngineR_firstzero_ac", (DL_FUNC) &_CompEngineR_firstzero_ac, 2},
    {"_CompEngineR_fluctanal_prop_r1", (DL_FUNC) &_CompEngineR_fluctanal_prop_r1, 1},
    {"_CompEngineR_histogram_mode", (DL_FUNC) &_CompEngineR_histogram_mode, 2},
    {"_CompEngineR_localsimple_taures", (DL_FUNC) &_CompEngineR_localsimple_taures, 2},
    {"_CompEngineR_motiftwo_entro3", (DL_FUNC) &_CompEngineR_motiftwo_entro3, 1},
    {"_CompEngineR_outlierinclude_mdrmd", (DL_FUNC) &_CompEngineR_outlierinclude_mdrmd, 2},
    {"_CompEngineR_sampen_first", (DL_FUNC) &_CompEngineR_sampen_first, 1},
    {"_CompEngineR_sampenc", (DL_FUNC) &_CompEngineR_sampenc, 3},
    {"_CompEngineR_spreadrandomlocal_meantaul", (DL_FUNC) &_CompEngineR_spreadrandomlocal_meantaul, 2},
    {"_CompEngineR_std1st_der", (DL_FUNC) &_CompEngineR_std1st_der, 1},
    {"_CompEngineR_trev_num", (DL_FUNC) &_CompEngineR_trev_num, 1},
    {"_CompEngineR_walker_propcross", (DL_FUNC) &_CompEngineR_walker_propcross, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_CompEngineR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
