#include <Rcpp.h>
#include "Ckmeans.1d.dp.h"
#include <vector>   
#include "R.h"
#include "Rmath.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector getPIP_kmeans(IntegerVector seeds, IntegerVector scans, NumericVector mzs, NumericVector ints, IntegerVector clu, double mztol, int gap, int min_width, int max_width, double alpha) {
  int idx = seeds.size();
  int end = mzs.size();
  int cid = 1;
  
  for (int i=0; i<idx; i++){
    int ref = seeds[i]-1;
    if (clu[ref]!=0){
      continue;
    }
    
    double refMz = mzs[ref];
    int refScan = scans[ref];
    double refInt = ints[ref];
    
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
      if (scans[j]>=leftScan && scans[j]<=rightScan){
        roi.push_back(j);
      }
    }
    for (int k=ref+1; k<end; k++){
      if (mzs[k]>downMz){
        break;
      }
      if (scans[k]>=leftScan && scans[k]<=rightScan){
        roi.push_back(k);
      }
    }
    
    NumericVector roi_mzs = mzs[roi];
    IntegerVector roi_scans = scans[roi];
    NumericVector roi_ints = ints[roi];
    
    // kmeans cluster
    NumericVector mzdiffs = (roi_mzs - refMz) * (roi_mzs - refMz);
    int roi_length = roi.size();
    std::vector<double> mzdiffs1 = Rcpp::as<std::vector<double> > (mzdiffs);
    std::vector<double> input(roi_length+1);
    std::vector<double> cluster(roi_length);
    
    for (size_t j=1; j<input.size(); ++j) {
      input[j] = mzdiffs1[j-1];
    }
    
    ClusterResult 
      result = kmeans_1d_dp(input, 1, 5);
    
    for (int k=1; k <= roi_length; ++k) {
      cluster[k-1] = (int) result.cluster[k];
    }
    
    IntegerVector sels;
    for (int l=0; l<cluster.size(); l++){
      if (cluster[l]==1){
        sels.push_back(l);
      }
    }
    
    // collect pic
    IntegerVector sel_id = roi[sels];
    NumericVector sel_mz = roi_mzs[sels];
    IntegerVector sel_scan = roi_scans[sels];
    NumericVector sel_ints = roi_ints[sels];
    
    int pic_leftscan = min(sel_scan);
    int pic_rightscan = max(sel_scan);
    
    if (pic_rightscan - pic_leftscan < min_width){
      clu[sel_id] = -1;
      continue;
    }
    
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
    
    for (int ss=refScan; ss>=pic_leftscan; ss--){
      if (gap0 > gap){
        break;
      }
      winner_int = -100000;
      winner = -1;
      w = alpha;
      at = 2*S1-S2;
      bt = w/(1-w)*(S1-S2);
      xhat = at+bt>0?at+bt:0;
      
      for (int j=0; j<sel_id.size(); j++){
        if (sel_scan[j]==ss && abs(sel_ints[j]-xhat)<abs(winner_int-xhat)){
          winner_int = sel_ints[j];
          winner = sel_id[j];
        }
      }
      
      if (winner>0){
        clu[winner] = cid;
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
      winner_int = -100000;
      winner = -1;
      w = alpha;
      at = 2*S1-S2;
      bt = w/(1-w)*(S1-S2);
      xhat = at+bt>0?at+bt:0;;
      
      for (int j=0; j<sel_id.size(); j++){
        if (sel_scan[j]==ss && abs(sel_ints[j]-xhat)<abs(winner_int-xhat)){
          winner_int = sel_ints[j];
          winner = sel_id[j];
        }
      }
      
      if (winner>0){
        clu[winner] = cid;
        xx = winner_int;
        gap0 = 0;
      } else {
        xx = xhat;
        gap0 = gap0 + 1;
      }
      S1 = w*xx+(1-w)*S1;
      S2 = w*S1+(1-w)*S2;
    }
    
    cid = cid + 1;
  }
  return (clu);
}