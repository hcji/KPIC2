#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector getPIP(IntegerVector seeds, IntegerVector scans, NumericVector mzs, IntegerVector clu, double mztol, int gap) {
  int idx = seeds.size();
  int end = mzs.size();
  int cid = 1;
  for (int i=0; i<idx; i++) {
    int ref = seeds[i]-1;
    if (clu[ref]!=0){
      continue;
    }
    
    clu[ref] = cid;
    double refMz = mzs[ref];
    int refScan = scans[ref];
    
    // left side
    int gap0 = 0;
    int np =1;
    int thisScan = refScan-1;
    double mu = refMz;
    double sigma = mztol/3;
    
    int uphit;
    int downhit;
    double updiff;
    double downdiff;
    
    while (true) {
      if (gap0>=gap) {
        break;
      }
      // find closest
      uphit = -1;
      downhit = -1;
      updiff = 10^5;
      downdiff = 10^5;
      
      for (int j = ref; j>0; j--){
        updiff = mu - mzs[j];
        if (updiff>3*sigma){
          break;
        }
        if (scans[j]==thisScan && clu[j]==0){
          uphit = j;
          break;
        }
      }
      
      
      for (int j = ref; j<end; j++){
        downdiff = mzs[j]-mu;
        if (downdiff>3*sigma){
          break;
        }
        if (scans[j]==thisScan && clu[j]==0){
          downhit = j;
          break;
        }
      }
      
      if (uphit<0 && downhit<0){
        gap0 = gap0 + 1;
        thisScan = thisScan - 1;
        continue;
      }
      
      // compare
      double wdiff = 0;
      int winner = -1;
      if (updiff <= downdiff){
        wdiff = updiff;
        winner = uphit;
      } else {
        wdiff = downdiff;
        winner = downhit;
      }
      
      // reset mu and sigma
      double this_mz = mzs[winner];
      mu = (np * mu + this_mz) / (np + 1);
      sigma = sqrt((np * sigma * sigma + wdiff * wdiff) / (np + 1));
      
      clu[winner] = cid;
      gap0 = 0;
      np = np + 1;
      thisScan = thisScan - 1;
    }
    
    
    // right side
    gap0 = 0;
    np = 1;
    thisScan = refScan+1;
    mu = refMz;
    sigma = mztol/3;
    while (true) {
      if (gap0>=gap) {
        break;
      }
      // find closest
      uphit = -1;
      downhit = -1;
      updiff = 10^5;
      downdiff = 10^5;
      
      for (int j = ref; j>0; j--){
        updiff = mu - mzs[j];
        if (updiff>3*sigma){
          break;
        }
        if (scans[j]==thisScan && clu[j]==0){
          uphit = j;
          break;
        }
      }
      
      
      for (int j = ref; j<end; j++){
        downdiff = mzs[j]-mu;
        if (downdiff>3*sigma){
          break;
        }
        if (scans[j]==thisScan && clu[j]==0){
          downhit = j;
          break;
        }
      }
      
      if (uphit<0 && downhit<0){
        gap0 = gap0 + 1;
        thisScan = thisScan + 1;
        continue;
      }
      
      // compare
      double wdiff = 0;
      int winner = -1;
      if (updiff <= downdiff){
        wdiff = updiff;
        winner = uphit;
      } else {
        wdiff = downdiff;
        winner = downhit;
      }
      
      // reset mu and sigma
      double this_mz = mzs[winner];
      mu = (np * mu + this_mz) / (np + 1);
      sigma = sqrt((np * sigma * sigma + wdiff * wdiff) / (np + 1));
      
      clu[winner] = cid;
      gap0 = 0;
      np = np + 1;
      thisScan = thisScan + 1;
    }
    
    cid = cid + 1;
  }
  return clu;
}

