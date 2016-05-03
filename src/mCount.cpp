#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List mCount(SEXP x, SEXP y, SEXP z, SEXP a, SEXP b) {
	
	Rcpp::NumericMatrix matrix1(x);
	Rcpp::NumericMatrix matrix2(y);
	Rcpp::NumericMatrix mscorematrix(z);
	int above = Rcpp::as<int>(a);
	int between = Rcpp::as<int>(b);
	int cols1 = matrix1.ncol();
	int rows1 = matrix1.nrow();
	int cols2 = matrix2.ncol();
	Rcpp::NumericMatrix count_matrix(rows1, cols2);
	Rcpp::NumericMatrix count_mscore_matrix(rows1, cols2);
	Rcpp::NumericVector temp_row1(cols1);
	Rcpp::NumericVector temp_row2(cols2);
	Rcpp::NumericVector mscore(rows1); 
	Rcpp::NumericVector morder(rows1);
	Rcpp::NumericVector mcount(rows1); 
	Rcpp::NumericVector mcutoff(rows1);
	Rcpp::NumericVector means2(rows1);
	int count, idx;
	double cutoff, sum2;
	List ret;
	
	//computing M^(j)_(group1) for each feature i
	for (int i = 0; i < rows1; i++) {
		temp_row1 = matrix1(i,_);
		temp_row2 = matrix2(i,_);
		std::sort(temp_row1.begin(), temp_row1.end(), std::greater<double>());
		std::sort(temp_row2.begin(), temp_row2.end(), std::greater<double>());
		matrix1(i,_) = temp_row1;
		matrix2(i,_) = temp_row2;
		
		sum2 = 0;
		for (int j = 0; j < cols2; j++) {
			count = 0;
			sum2 = sum2 + matrix2(i,j);
			cutoff = matrix2(i,j) + between - 1;
			for (int k = 0; k < cols1; k++) {
				if (matrix1(i,k) >= cutoff && matrix1(i,k) >= above) {
					count++;
				} else {
					break;
				}
			}
			count_matrix(i,j) = count;
			count_mscore_matrix(i,j) = mscorematrix(j,count);
		}
		means2(i) = sum2 / cols2;	
		idx = which_min(count_mscore_matrix(i,_));
		mscore(i) = count_mscore_matrix(i,idx);
		morder(i) = idx;
		mcount(i) = count_matrix(i,idx);
		mcutoff(i) = matrix2(i,idx) + between - 1;
	}
	
	ret["mMs"] = mscore;
	ret["order"] = morder;
	ret["count"] = mcount;
	ret["cutoff"] = mcutoff;
	ret["means2"] = means2;
	return ret;
}