// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cerf
double cerf(double x);
RcppExport SEXP ProFit_cerf(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    __result = Rcpp::wrap(cerf(x));
    return __result;
END_RCPP
}
// profitMakeSersic
NumericMatrix profitMakeSersic(const IntegerMatrix& CALCREGION, const double XCEN, const double YCEN, const double MAG, const double RE, const double NSER, const double ANG, const double AXRAT, const double BOX, const double MAGZERO, const bool ROUGH, const NumericVector& XLIM, const NumericVector& YLIM, const IntegerVector& DIM, const int UPSCALE, const int MAXDEPTH, const double RESWITCH, const double ACC, const bool DOCALCREGION, const double REMAX);
RcppExport SEXP ProFit_profitMakeSersic(SEXP CALCREGIONSEXP, SEXP XCENSEXP, SEXP YCENSEXP, SEXP MAGSEXP, SEXP RESEXP, SEXP NSERSEXP, SEXP ANGSEXP, SEXP AXRATSEXP, SEXP BOXSEXP, SEXP MAGZEROSEXP, SEXP ROUGHSEXP, SEXP XLIMSEXP, SEXP YLIMSEXP, SEXP DIMSEXP, SEXP UPSCALESEXP, SEXP MAXDEPTHSEXP, SEXP RESWITCHSEXP, SEXP ACCSEXP, SEXP DOCALCREGIONSEXP, SEXP REMAXSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type CALCREGION(CALCREGIONSEXP);
    Rcpp::traits::input_parameter< const double >::type XCEN(XCENSEXP);
    Rcpp::traits::input_parameter< const double >::type YCEN(YCENSEXP);
    Rcpp::traits::input_parameter< const double >::type MAG(MAGSEXP);
    Rcpp::traits::input_parameter< const double >::type RE(RESEXP);
    Rcpp::traits::input_parameter< const double >::type NSER(NSERSEXP);
    Rcpp::traits::input_parameter< const double >::type ANG(ANGSEXP);
    Rcpp::traits::input_parameter< const double >::type AXRAT(AXRATSEXP);
    Rcpp::traits::input_parameter< const double >::type BOX(BOXSEXP);
    Rcpp::traits::input_parameter< const double >::type MAGZERO(MAGZEROSEXP);
    Rcpp::traits::input_parameter< const bool >::type ROUGH(ROUGHSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type XLIM(XLIMSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type YLIM(YLIMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type DIM(DIMSEXP);
    Rcpp::traits::input_parameter< const int >::type UPSCALE(UPSCALESEXP);
    Rcpp::traits::input_parameter< const int >::type MAXDEPTH(MAXDEPTHSEXP);
    Rcpp::traits::input_parameter< const double >::type RESWITCH(RESWITCHSEXP);
    Rcpp::traits::input_parameter< const double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< const bool >::type DOCALCREGION(DOCALCREGIONSEXP);
    Rcpp::traits::input_parameter< const double >::type REMAX(REMAXSEXP);
    __result = Rcpp::wrap(profitMakeSersic(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION, REMAX));
    return __result;
END_RCPP
}
// profitBruteConv
NumericMatrix profitBruteConv(const NumericMatrix& IMG, const NumericMatrix& PSF, const IntegerMatrix& CALCREGION, const bool DOCALCREGION);
RcppExport SEXP ProFit_profitBruteConv(SEXP IMGSEXP, SEXP PSFSEXP, SEXP CALCREGIONSEXP, SEXP DOCALCREGIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type IMG(IMGSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type PSF(PSFSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type CALCREGION(CALCREGIONSEXP);
    Rcpp::traits::input_parameter< const bool >::type DOCALCREGION(DOCALCREGIONSEXP);
    __result = Rcpp::wrap(profitBruteConv(IMG, PSF, CALCREGION, DOCALCREGION));
    return __result;
END_RCPP
}
// profitBruteConv2
NumericMatrix profitBruteConv2(const NumericMatrix& IMG, const NumericMatrix& PSF, const IntegerMatrix& CALCREGION, const bool DOCALCREGION);
RcppExport SEXP ProFit_profitBruteConv2(SEXP IMGSEXP, SEXP PSFSEXP, SEXP CALCREGIONSEXP, SEXP DOCALCREGIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type IMG(IMGSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type PSF(PSFSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type CALCREGION(CALCREGIONSEXP);
    Rcpp::traits::input_parameter< const bool >::type DOCALCREGION(DOCALCREGIONSEXP);
    __result = Rcpp::wrap(profitBruteConv2(IMG, PSF, CALCREGION, DOCALCREGION));
    return __result;
END_RCPP
}
// profitDownsample
NumericMatrix profitDownsample(const NumericMatrix& IMG, const int DOWNSAMPLEFAC);
RcppExport SEXP ProFit_profitDownsample(SEXP IMGSEXP, SEXP DOWNSAMPLEFACSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type IMG(IMGSEXP);
    Rcpp::traits::input_parameter< const int >::type DOWNSAMPLEFAC(DOWNSAMPLEFACSEXP);
    __result = Rcpp::wrap(profitDownsample(IMG, DOWNSAMPLEFAC));
    return __result;
END_RCPP
}
// profitUpsample
NumericMatrix profitUpsample(const NumericMatrix& IMG, const int UPSAMPLEFAC);
RcppExport SEXP ProFit_profitUpsample(SEXP IMGSEXP, SEXP UPSAMPLEFACSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type IMG(IMGSEXP);
    Rcpp::traits::input_parameter< const int >::type UPSAMPLEFAC(UPSAMPLEFACSEXP);
    __result = Rcpp::wrap(profitUpsample(IMG, UPSAMPLEFAC));
    return __result;
END_RCPP
}
