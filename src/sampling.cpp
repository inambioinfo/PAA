#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List sampling(SEXP x, SEXP c, SEXP c1, SEXP c2, SEXP tc1, SEXP tc2) {
	
	Rcpp::NumericMatrix datamatrix(x);
	Rcpp::CharacterVector cols(c);
	Rcpp::CharacterVector cols1(c1);
	Rcpp::CharacterVector cols2(c2);
	Rcpp::IntegerVector tcols1(tc1);
	Rcpp::IntegerVector tcols2(tc2);
	
	int p = datamatrix.nrow();
	int n1 = cols1.size();
	int n2 = cols2.size();
	int n = cols.size();
	int tn1 = tcols1.size();
	int tn2 = tcols2.size();
	int tn = tn1 + tn2;
	int idx_tcols1, idx_tcols2, idx_traincols, idx_traincols1, idx_traincols2, idx_testcols, idx_testcols1, idx_testcols2;
	
	Rcpp::NumericMatrix trainmatrix(p,n-tn);
	Rcpp::NumericMatrix trainmatrix1(p,n1-tn1);
	Rcpp::NumericMatrix trainmatrix2(p,n2-tn2);
	Rcpp::NumericMatrix testmatrix(p,tn);
	Rcpp::CharacterVector traincols(n-tn);
	Rcpp::CharacterVector traincols1(n1-tn1);
	Rcpp::CharacterVector traincols2(n2-tn2);
	Rcpp::CharacterVector testcols(tn);
	Rcpp::CharacterVector testcols1(tn1);
	Rcpp::CharacterVector testcols2(tn2);
	List ret;
	
	idx_tcols1 = idx_tcols2 = idx_traincols = idx_traincols1 = idx_traincols2 = idx_testcols = idx_testcols1 = idx_testcols2 = 0;
	
	for (int i = 0; i < n; i++) {
		if (i == tcols1(idx_tcols1)) {
			testmatrix(_,idx_testcols) = datamatrix(_,i);
			testcols(idx_testcols) = cols(i);
			testcols1(idx_testcols1) = cols(i);
			if (idx_tcols1 < tn1-1) {
				idx_tcols1++;
				idx_testcols1++;
			}
			idx_testcols++;
		} else if (i == tcols2(idx_tcols2)) {
			testmatrix(_,idx_testcols) = datamatrix(_,i);
			testcols(idx_testcols) = cols(i);
			testcols2(idx_testcols2) = cols(i);
			if (idx_tcols2 < tn2-1) {
				idx_tcols2++;
				idx_testcols2++;
			}
			idx_testcols++;
		} else if (i < n1) {
			trainmatrix(_,idx_traincols) = datamatrix(_,i);
			trainmatrix1(_,idx_traincols1) = datamatrix(_,i);
			traincols(idx_traincols) = cols(i);
			traincols1(idx_traincols1) = cols(i);
			idx_traincols++;
			idx_traincols1++;
		} else {
			trainmatrix(_,idx_traincols) = datamatrix(_,i);
			trainmatrix2(_,idx_traincols2) = datamatrix(_,i);
			traincols(idx_traincols) = cols(i);
			traincols2(idx_traincols2) = cols(i);
			idx_traincols++;
			idx_traincols2++;
		}
	}
	
	ret["trainmatrix"] = trainmatrix;
	ret["trainmatrix1"] = trainmatrix1;
	ret["trainmatrix2"] = trainmatrix2;
	ret["testmatrix"] = testmatrix;
	ret["traincols"] = traincols;
	ret["traincols1"] = traincols1;
	ret["traincols2"] = traincols2;
	ret["testcols"] = testcols;
	ret["testcols1"] = testcols1;
	ret["testcols2"] = testcols2;
	return ret;
}