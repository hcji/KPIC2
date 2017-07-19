//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix waveft(NumericVector omega,NumericVector scales) {
  double StpFrq = omega[1];
  int NbFrq  = omega.size();
  double SqrtNbFrq = sqrt(NbFrq);
  double cfsNORM = sqrt(StpFrq)*SqrtNbFrq;
  int NbSc = scales.size();
  NumericMatrix wft(NbSc,NbFrq);
  NumericVector mul = sqrt(scales/1.32934)*cfsNORM;

  for (int jj=0; jj<NbSc; jj++){
    NumericVector scapowered = (scales[jj]*omega);
    NumericVector expnt = -pow(scapowered,2)/2;
    wft.row(jj) = mul[jj]*pow(scapowered,2)*exp(expnt);
  }
  return(wft);
}

// [[Rcpp::export]]
List cwtft(NumericVector val) {
  double meanSIG = mean(val);
  vec xx = val-meanSIG;
  int n = xx.size();
  double dt = 2;

  NumericVector omega;
  for (int i=floor(n/2); i>0; i--){
    omega.push_front(i);
  }
  for (int i=floor((n-1)/2); i>0; i--){
    omega.push_back(-i);
  }
  omega.push_front(0);
  omega = omega*((2*3.1415926)/(n*dt));

  cx_vec ff = fft(xx);

  double s0 = dt;
  double ds = 0.2;
  int NbSc = floor(log2(n)/ds);

  vec scales1 = linspace<vec> (0,(NbSc-1),NbSc);
  scales1 = s0*exp2((scales1)*ds);
  NumericVector scales = wrap(scales1);
  scales.push_front(1);

  NumericMatrix psift = waveft(omega,scales);

  cx_mat MAT(NbSc+1, ff.size());
  for (int i=0; i<=NbSc; i++) {
    MAT.row(i) = ff.t();
  }

  cx_mat cwtcfs = MAT % as<mat>(psift);
  cwtcfs = cwtcfs.t();
  cwtcfs = ifft(cwtcfs)/ff.size();
  mat cwtcfs1 = real(cwtcfs.t());

  List output;
  output["scales"] = scales;
  output["cwt2d"] = cwtcfs1;

  return output;
}

// [[Rcpp::export]]
LogicalMatrix localMax(NumericMatrix cwt2d) {
  int cols = cwt2d.ncol();
  int rows = cwt2d.nrow();
  LogicalMatrix result(rows,cols);
  NumericVector rowi;
  for (int i=0; i<rows; i++){
    rowi = cwt2d.row(i);
    rowi.push_back(0);rowi.push_back(0);
    rowi.push_front(0);rowi.push_front(0);
    for (int j=2; j<cols-2; j++){
      if (rowi[j]>rowi[j+1] && rowi[j]>rowi[j+2] && rowi[j]>rowi[j-1] && rowi[j]>rowi[j-2]){
        result(i,j-2) = true;
      } else {
        result(i,j-2) = false;
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
LogicalMatrix localMin(NumericMatrix cwt2d) {
  int cols = cwt2d.ncol();
  int rows = cwt2d.nrow();
  LogicalMatrix result(rows,cols);
  NumericVector rowi;
  for (int i=0; i<rows; i++){
    rowi = cwt2d.row(i);
    rowi.push_back(0);rowi.push_back(0);
    rowi.push_front(0);rowi.push_front(0);
    for (int j=2; j<cols-2; j++){
      if (rowi[j]<rowi[j+1] && rowi[j]<rowi[j+2] && rowi[j]<rowi[j-1] && rowi[j]<rowi[j-2]){
        result(i,j-2) = true;
      } else {
        result(i,j-2) = false;
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
List ridgesDetection(NumericMatrix cwt2d, NumericVector val) {
  List output;
  int n_rows = cwt2d.nrow();
  int n_cols = cwt2d.ncol();
  LogicalMatrix local_max = localMax(cwt2d);

  NumericVector cols_small_peaks;
  NumericVector rows_small_peaks;

  for (int j=0; j<n_cols; j++){
    for (int i=0; i<5; i++){
      if (local_max(i,j)){
        rows_small_peaks.push_back(i);
        cols_small_peaks.push_back(j);
        break;
      }
    }
  }

  List ridges_rows;
  List ridges_cols;
  for (int i=0; i<cols_small_peaks.size(); i++){
    NumericVector rows;
    NumericVector cols;
    cols.push_back(cols_small_peaks(i));
    rows.push_back(rows_small_peaks(i));

    int row_best = rows_small_peaks(i);
    int col_plus = cols_small_peaks(i);
    int col_minus = cols_small_peaks(i);

    for (int j=1; j<n_rows; j++){
      int row_plus = row_best + j;
      int row_minus = row_best - j;

      if (row_minus>0 && 2<=col_minus && col_minus<=n_cols-2){
        if (local_max(row_minus, col_minus + 1)){col_minus = col_minus + 1;
        }else if(local_max(row_minus, col_minus + 2)){col_minus = col_minus + 2;
        }else if(local_max(row_minus, col_minus - 2)){col_minus = col_minus - 2;
        }else if(local_max(row_minus, col_minus - 1)){col_minus = col_minus - 1;
        }else if(local_max(row_minus, col_minus)){col_minus = col_minus;
        }else{col_minus = -1;}
        if (col_minus != -1){
          rows.push_front(row_minus);
          cols.push_front(col_minus);
        }
      }

      if (row_plus<n_rows && 2<=col_plus && col_plus<=n_cols-2){
        if (local_max(row_plus, col_plus + 1)){col_plus = col_plus + 1;
        }else if(local_max(row_plus, col_plus + 2)){col_plus = col_plus + 2;
        }else if(local_max(row_plus, col_plus - 2)){col_plus = col_plus - 2;
        }else if(local_max(row_plus, col_plus - 1)){col_plus = col_plus - 1;
        }else if(local_max(row_plus, col_plus)){col_plus = col_plus;
        }else{col_plus = -1;}
        if (col_plus != -1){
          rows.push_back(row_plus);
          cols.push_back(col_plus);
        }
      }
    }

    if (rows.size()>=7){
      ridges_rows.push_back(rows);
      ridges_cols.push_back(cols);
    }
  }

  output["ridges_rows"] = ridges_rows;
  output["ridges_cols"] = ridges_cols;

  return output;
}

// [[Rcpp::export]]
NumericVector peaksPosition(NumericVector val, List ridges, NumericMatrix cwt2d){
  List output;
  List ridges_rows = ridges["ridges_rows"];
  List ridges_cols = ridges["ridges_cols"];
  List ridges_rows_select;
  List ridges_cols_select;
  NumericVector peaks;

  int n_cols = cwt2d.ncol();
  int n_rows = cwt2d.nrow();

  LogicalMatrix local_minus = localMin(cwt2d);
  LogicalMatrix neg(n_rows,n_cols);
  for (int i=0; i<n_rows; i++){
    LogicalVector negi = cwt2d.row(i)<0;
    negi = negi|local_minus.row(i);
    neg.row(i) = negi;
  }
  neg.column(0) <- true;
  neg.column(neg.ncol()-1) <- true;

  for (int i=0; i<ridges_rows.size(); i++){
    rowvec ridge_rows = as<rowvec>(ridges_rows(i));
    rowvec ridge_cols = as<rowvec>(ridges_cols(i));

    NumericVector inds;
    for (int j=0; j<ridge_rows.size(); j++){
      if (cwt2d(ridge_rows(j), ridge_cols(j)) > 0){
        inds.push_back(j);
      }
    }

    if (inds.size()>0){
      int col = round(mean(ridge_cols.elem(as<uvec>(inds))));
      int row = ridge_rows(abs(ridge_cols-col).index_min());
      int peak = col;
      double top = 0;

      for (int j=0; j<n_cols; j++){
        int col_minus = col-j;
        if (col_minus >= 0){
          if (val(col_minus) > top){
            peak = col_minus;
            top = val(col_minus);
          }
          if (neg(row,col_minus)){
            break;}
        } else {break;}
      }
      for (int j=0; j<n_cols; j++){
        int col_plus = col+j;
        if (col_plus < n_cols){
          if (val(col_plus) > top){
            peak = col_plus;
            top = val(col_plus);
          }
          if (neg(row, col_plus)){
            break;}
        } else {break;}
      }
      peaks.push_back(peak);
    } else {
      int l_col = ridge_cols.size()-1;
      int cols_start = n_cols;
      int cols_end = 0;
      for (int j=0; j<l_col; j++){
        if (ridge_cols(j) > cols_end) {cols_end = ridge_cols(j);}
        if (ridge_cols(j) < cols_start) {cols_start = ridge_cols(j);}
      }
      if (cols_start - 3 < 0){
        cols_start = 0;
      } else {cols_start = cols_start - 3;}
      if (cols_end + 4 > n_cols - 1){
        cols_end = n_cols - 1;
      } else {cols_end = cols_end + 4;}

      int peak = 0;
      double top = 0;
      for (int j=cols_start; j<cols_end; j++){
        if (val(j) >= top){
          peak = j;
          top = val(j);}
      }
      peaks.push_back(peak);
    }
  }
  return peaks;
}

// [[Rcpp::export]]
List getSignal(NumericMatrix cwt2d, List ridges, NumericVector peaks){
  List output;
  NumericVector signals;
  NumericVector ridge_len;
  List ridges_rows = ridges["ridges_rows"];
  List ridges_cols = ridges["ridges_cols"];

  for (int ind=0; ind<peaks.size(); ind++){
    NumericVector ridge_rows = ridges_rows(ind);
    NumericVector ridge_cols = ridges_cols(ind);

    double signal = 0;
    int best_len;
    for (int i=0; i<ridge_rows.size(); i++){
      int a = ridge_rows(i);
      int b = ridge_cols(i);
      double s = cwt2d(a, b);
      if (s >= signal){
        signal = s;
        best_len = a;}
    }
    ridge_len.push_back(best_len);
    signals.push_back(signal);
  }

  output["signals"] = signals;
  output["ridge_lens"] = ridge_len;
  return output;
}

