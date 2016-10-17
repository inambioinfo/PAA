#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::NumericMatrix mMsMatrix(SEXP x, SEXP y) {
	
	int n1 = Rcpp::as<int>(x);
	int n2 = Rcpp::as<int>(y);
	Rcpp::NumericMatrix P_Mij_matrix2(n1, n2+1);
	Rcpp::NumericMatrix pij_matrix(n1, n2+1);
	double p;
	
	for (int i = 0; i < n1; i++) {
		for (int m = 0; m < n2+1; m++) { 
			P_Mij_matrix2(i,m) = R::choose(n1+n2-m-i-1, n2-i-1)*R::choose(m+i,i)/R::choose(n1+n2, n1);
			//choose(n1+n2-m-i, n2-i)*choose(m+i-1,i-1)/choose(n1+n2, n1) -> because of C: i = (i+1)
		}
	}
	
	for (int i = 0; i < n1; i++) {
		for (int m = 0; m < n2+1; m++) {
			p = 0;
			for (int k = m; k < n2+1; k++) {
				p = p + P_Mij_matrix2(i,k);
			}
			pij_matrix(i,m) = p;
		}
	}
	
	return pij_matrix;
}