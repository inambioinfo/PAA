// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// joinMCountResults
List joinMCountResults(SEXP p1, SEXP p2, SEXP m1, SEXP m2, SEXP l1, SEXP l2);
RcppExport SEXP PAA_joinMCountResults(SEXP p1SEXP, SEXP p2SEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP l1SEXP, SEXP l2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type p1(p1SEXP );
        Rcpp::traits::input_parameter< SEXP >::type p2(p2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type m1(m1SEXP );
        Rcpp::traits::input_parameter< SEXP >::type m2(m2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type l1(l1SEXP );
        Rcpp::traits::input_parameter< SEXP >::type l2(l2SEXP );
        List __result = joinMCountResults(p1, p2, m1, m2, l1, l2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mCount
List mCount(SEXP x, SEXP y, SEXP z, SEXP a, SEXP b);
RcppExport SEXP PAA_mCount(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type x(xSEXP );
        Rcpp::traits::input_parameter< SEXP >::type y(ySEXP );
        Rcpp::traits::input_parameter< SEXP >::type z(zSEXP );
        Rcpp::traits::input_parameter< SEXP >::type a(aSEXP );
        Rcpp::traits::input_parameter< SEXP >::type b(bSEXP );
        List __result = mCount(x, y, z, a, b);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mMsMatrix
Rcpp::NumericMatrix mMsMatrix(SEXP x, SEXP y);
RcppExport SEXP PAA_mMsMatrix(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type x(xSEXP );
        Rcpp::traits::input_parameter< SEXP >::type y(ySEXP );
        Rcpp::NumericMatrix __result = mMsMatrix(x, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sampling
List sampling(SEXP x, SEXP c, SEXP c1, SEXP c2, SEXP tc1, SEXP tc2);
RcppExport SEXP PAA_sampling(SEXP xSEXP, SEXP cSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP tc1SEXP, SEXP tc2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type x(xSEXP );
        Rcpp::traits::input_parameter< SEXP >::type c(cSEXP );
        Rcpp::traits::input_parameter< SEXP >::type c1(c1SEXP );
        Rcpp::traits::input_parameter< SEXP >::type c2(c2SEXP );
        Rcpp::traits::input_parameter< SEXP >::type tc1(tc1SEXP );
        Rcpp::traits::input_parameter< SEXP >::type tc2(tc2SEXP );
        List __result = sampling(x, c, c1, c2, tc1, tc2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}