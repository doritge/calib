#include <Rcpp.h>
using namespace Rcpp;

int check_index(int n, int i)
  // wrap index if beyond limits
{
  if (i < 0)
    return n + i;
  if(i >= n)
    return i % n;
  return i;
}

// [[Rcpp::export]]
NumericMatrix iterate_cyclic(NumericMatrix X,
                             Rcpp::DataFrame P,
                             int states,
                             int threshold){
  int m = X.nrow();
  int n = X.ncol();
  int k = P.nrows();

  IntegerVector dx = P["x"];
  IntegerVector dy = P["y"];
  IntegerVector w = P["weight"];

  NumericMatrix X_new = clone(X);

  for(int y = 0; y < m; y++) {
    for(int x = 0; x < n; x++){
      int count = 0;
      int v = X(y,x);
      int v_n = (v+1) % states;
      for(int z = 0; z < k; z++){
        int iy  = check_index(m, y + dx[z]);
        int ix  = check_index(n, x + dy[z]);
        if (X(iy, ix) * w[z] == v_n){
          count++;
        }
      }
      if (count >= threshold){
        X_new(y,x) = v_n;
      }
    }
  }
  return X_new;
}

// [[Rcpp::export]]
IntegerMatrix iterate_index(IntegerMatrix X,
                             Rcpp::DataFrame P,
                             IntegerVector R){
  int m = X.nrow();
  int n = X.ncol();
  //int states = indices.size();

  IntegerVector dx = P["x"];
  IntegerVector dy = P["y"];
  IntegerVector w = P["weight"];

  IntegerMatrix X_new = clone(X);

  for(int y = 0; y < m; y++) {
    for(int x = 0; x < n; x++){

      int index = R[X(y,x)];
      int iy  = check_index(m, y + dy[index]);
      int ix  = check_index(n, x + dx[index]);
      X_new(y, x) = X(iy, ix);
    }
  }
  return X_new;
}

// [[Rcpp::export]]
NumericMatrix iterate_life(NumericMatrix X,
                           Rcpp::DataFrame P,
                           int rule){
  // Game of life rule B3S23. 2 states Moore neighborhood
  int m = X.nrow();
  int n = X.ncol();
  int k = P.nrows();

  IntegerVector dx = P["x"];
  IntegerVector dy = P["y"];

  NumericMatrix X_new = clone(X);

  for(int y = 0; y < m; y++) {
    for(int x = 0; x < n; x++){
      int v = X(y, x);
      int count = 0;
      for(int z = 0; z < k; z++){
        int iy  = check_index(m, y + dx[z]);
        int ix  = check_index(n, x + dy[z]);
        count = count + X(iy, ix);
      }
      switch(rule){
      case 1:    // B3S23
        if (v == 0 && count == 3) X_new(y, x) = 1;   // birth
        if (v == 1 && !(count == 2 || count == 3)) X_new(y, x) = 0;    // death
      }
    }
  }
  return X_new;
}

// [[Rcpp::export]]
NumericMatrix iterate_total(NumericMatrix X,
                            Rcpp::DataFrame P,
                            IntegerVector R){
  int m = X.nrow();
  int n = X.ncol();
  int k = P.nrows();

  IntegerVector dx = P["x"];
  IntegerVector dy = P["y"];
  IntegerVector w = P["weight"];

  NumericMatrix X_new = clone(X);

  for(int y = 0; y < m; y++) {
    for(int x = 0; x < n; x++){
      int count = 0;
      for(int z = 0; z < k; z++){
        int iy  = check_index(m, y + dx[z]);
        int ix  = check_index(n, x + dy[z]);
        count = count + X(iy, ix) * w[z];
      }
      X_new(y, x) = R[count];
    }
  }
  return X_new;
}
