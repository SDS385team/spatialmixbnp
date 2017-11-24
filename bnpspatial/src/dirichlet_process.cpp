#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector unifsb_dp(
    const size_t n, 
    const double concentration) {
  // Contains result
  NumericVector sbdp(n);
  // Draws v from beta distribution and 1 - v
  NumericVector v = rbeta(n, 1, concentration);
  NumericVector prod1v = cumprod(1-v);
  // Stick breaking
  sbdp[0] = v[0];
  for (size_t i = 1; i < n; i++) {
    sbdp[i] = prod1v[i] * v[i];
  }
  return sbdp;
}
  
// [[Rcpp::export]]
NumericVector polya_dp(
    const size_t n, 
    const double concentration) {
  // Contains result
  NumericVector sbdp(n);
  // Draws v from beta distribution and 1 - v
  NumericVector v = rbeta(n, 1, concentration);
  NumericVector prod1v = cumprod(1-v);
  // Stick breaking
  sbdp[0] = v[0];
  for (size_t i = 1; i < n; i++) {
    sbdp[i] = prod1v[i] * v[i];
  }
  return sbdp;
}
