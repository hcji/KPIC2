#include <Rcpp.h>
#include <vector>
#include "R.h"
#include "Rmath.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector getROI(int seed, IntegerVector scans, NumericVector mzs, NumericVector ints, LogicalVector notused, double mztol, int max_width) {
  int ref = seed-1;
  int end = mzs.size();
  double refMz = mzs[ref];
  int refScan = scans[ref];

  // locate roi
  int leftScan = refScan - 0.5 * max_width;
  int rightScan = refScan + 0.5 * max_width;
  double upMz = refMz - mztol;
  double downMz = refMz + mztol;

  IntegerVector roi = Rcpp::IntegerVector::create(ref);
  for (int j=ref-1; j>0; j--){
    if (mzs[j]<upMz){
      break;
    }
    if (scans[j]>=leftScan && scans[j]<=rightScan && notused[j]){
      roi.push_back(j);
    }
  }
  for (int k=ref+1; k<end; k++){
    if (mzs[k]>downMz){
      break;
    }
    if (scans[k]>=leftScan && scans[k]<=rightScan && notused[k]){
      roi.push_back(k);
    }
  }

  return roi;
}

// [[Rcpp::export]]
IntegerVector collectPIC(int refScan, double refMz, double refInt, IntegerVector sel_id, IntegerVector sel_scan, NumericVector sel_mz, NumericVector sel_ints, int gap, double alpha) {
  int pic_leftscan = min(sel_scan);
  int pic_rightscan = max(sel_scan);

  double xx = refInt;
  double xhat = refInt;
  double S1 = refInt;
  double S2 = refInt;
  double w;
  double at;
  double bt;
  int winner;
  double winner_int;
  int gap0 = 0;

  IntegerVector pic_id;

  for (int ss=refScan; ss>=pic_leftscan; ss--){
    if (gap0 > gap){
      break;
    }
    winner_int = -DOUBLE_XMAX;
    winner = -1;
    w = alpha;
    at = 2*S1-S2;
    bt = w/(1-w)*(S1-S2);
    xhat = at+bt>0?at+bt:0;

    for (int j=0; j<sel_id.size(); j++){
      if (sel_scan[j]==ss && fabs(sel_ints[j]-xhat)<fabs(winner_int-xhat)){
        winner_int = sel_ints[j];
        winner = sel_id[j];
      }
    }

    if (winner>0){
      pic_id.push_back(winner);
      xx = winner_int;
      gap0 = 0;
    } else {
      xx = xhat;
      gap0 = gap0 + 1;
    }
    S1 = w*xx+(1-w)*S1;
    S2 = w*S1+(1-w)*S2;
  }

  xx = refInt;
  xhat = refInt;
  S1 = refInt;
  S2 = refInt;
  gap0 = 0;
  for (int ss=refScan; ss<=pic_rightscan; ss++){
    if (gap0 > gap){
      break;
    }
    winner_int = -DOUBLE_XMAX;
    winner = -1;
    w = alpha;
    at = 2*S1-S2;
    bt = w/(1-w)*(S1-S2);
    xhat = at+bt>0?at+bt:0;

    for (int j=0; j<sel_id.size(); j++){
      if (sel_scan[j]==ss && fabs(sel_ints[j]-xhat)<fabs(winner_int-xhat)){
        winner_int = sel_ints[j];
        winner = sel_id[j];
      }
    }

    if (winner>0){
      pic_id.push_back(winner);
      xx = winner_int;
      gap0 = 0;
    } else {
      xx = xhat;
      gap0 = gap0 + 1;
    }
    S1 = w*xx+(1-w)*S1;
    S2 = w*S1+(1-w)*S2;
  }

  return pic_id;
}
