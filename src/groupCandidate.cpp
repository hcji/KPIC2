#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector findCandidate(int ind, double mzmin, double mzmax, double rtmin, double rtmax, NumericVector mz, NumericVector rt, NumericVector group) {
  NumericVector output;
  ind = ind - 1;
  output.push_back(ind);
  for (int i=1; i<mz.size(); i++){
    int lb = ind - i;
    if (lb >= 0){
      if (mz(lb) < mzmin) {break;}
      if (rt(lb) >= rtmin && rt(lb) <= rtmax && group(lb) == 0){
        output.push_back(lb);
      }
    } else {break;}
  }
  for (int i=1; i<mz.size(); i++){
    int rb = ind + i;
    if (rb < mz.size()){
      if (mz(rb) > mzmax) {break;}
      if (rt(rb) >= rtmin && rt(rb) <= rtmax && group(rb) == 0){
        output.push_back(rb);
      }
    } else {break;}
  }
  return output+1;
}
