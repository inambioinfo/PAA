#include <Rcpp.h>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List joinMCountResults(SEXP p1, SEXP p2, SEXP m1, SEXP m2, SEXP l1, SEXP l2) {
	
	Rcpp::NumericVector mscore1(p1);
	Rcpp::NumericVector mscore2(p2);
	Rcpp::NumericVector mean1(m1);
	Rcpp::NumericVector mean2(m2);
	string lab1 = Rcpp::as<std::string>(l1);
	string lab2 = Rcpp::as<std::string>(l2);
	
	int n = mean1.size();
	
	Rcpp::NumericVector p(n);
	Rcpp::NumericVector folds(n);
	Rcpp::CharacterVector labels(n);
	List ret;
	
	for (int i = 0; i < n; i++) {
		if (mscore1(i) <= mscore2(i)) {
			p(i) = mscore1(i);
			labels(i) = lab1;
		} else {
			p(i) = mscore2(i);
			labels(i) = lab2;
		}
		folds(i) = mean1(i) / mean2(i);
	}
	
	
	ret["p"] = p;
	ret["folds"] = folds;
	ret["labels"] = labels;
	return ret;
}