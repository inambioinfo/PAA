// This init file was created to register PAA native routines

#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP PAA_joinMCountResults(SEXP p1SEXP, SEXP p2SEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP l1SEXP, SEXP l2SEXP);
SEXP PAA_mCount(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP aSEXP, SEXP bSEXP);
SEXP PAA_mMsMatrix(SEXP xSEXP, SEXP ySEXP);
SEXP PAA_sampling(SEXP xSEXP, SEXP cSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP tc1SEXP, SEXP tc2SEXP);


R_CallMethodDef callMethods[] = {
    {"C_joinMCountResults", (DL_FUNC) &PAA_joinMCountResults, 6},
    {"C_mCount", (DL_FUNC) &PAA_mCount, 5},
    {"C_mMsMatrix", (DL_FUNC) &PAA_mMsMatrix, 2},
    {"C_sampling", (DL_FUNC) &PAA_sampling, 6},
    {NULL, NULL, 0}
};

void R_init_PAA(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
