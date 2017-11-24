#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES
#include <omp.h>
#include <RcppArmadillo.h>
#include <cmath>
#include <thread>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::plugins(openmp)]]

// inline void detect_and_set_cores() {
//   unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
//   if (concurentThreadsSupported != 0) {
//     Rcout << concurentThreadsSupported << endl;
//     omp_set_num_threads(concurentThreadsSupported);
//   }
// }

// This inline function creates one sample of multivariate t-distribution
// [[Rcpp::export]]
inline arma::vec rSt(
    const arma::vec& mu,
    const arma::mat& Sigma,
    const double df) {
  
  int d = mu.size();
  vec normal =  chol(Sigma).t() * randn(d);
  double chi2 = rchisq(1, df)[0];
  return mu + normal / sqrt(chi2 / df);
}

// This inline function evaluates the density of t distribution
// [[Rcpp::export]]
inline double dSt(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df) {
  
  int d = mu.n_elem;
  vec xcentered = x - mu;
  double innerterm =  - 0.5 * (df + d) * 
    log(1.0 + as_scalar(xcentered.t() * inv_sympd(Sigma) * xcentered) / df);
  double ldet;
  double sign;
  log_det(ldet, sign, Sigma); // compute and store the logarithm of determinant
  double extterm = log(tgamma(0.5 * (df + d))) - log(tgamma(0.5 * df)) - 0.5 * d * log(df * M_PI) - 0.5 * ldet;
  return exp(extterm + innerterm);
}


// [[Rcpp::export]]
int dp_gaussian_mixture(
    const arma::mat& y,
    const double alpha,
    const arma::vec& lambda,
    const arma::mat& Sigma,
    const double kappa,
    const double nu,
    const arma::mat& Omega
    ) {
  
  // Dimensions of the problem
  int T = y.n_rows;
  int d = y.n_cols;
  
  // Initialization of essential state particles
  uvec kt(T); // Allocation of observation at time t
  kt[0] = 0;
  int mt = 1; // Number of total components at time t, no need to store
  NumericVector ntj(T);
  ntj[0] = 1; // Counts of observations per component, T is the max number of clusters
  mat mutj(d, T); // Means per component
  mutj.col(0) = y.row(0).t();
  cube Stj(d, d, T); // Sum of squares per component
  Stj.slice(0) = mat(2, 2, fill::zeros);
  // Rcout << "kt: " << kt << endl << "mt: " << mt << endl << "ntj:" << ntj[0] <<
  //   endl << "mutj: " << mutj << endl << "ssj: " << ssj << endl;
  
  // Prior sampling constants
  double c0 = 2 * nu - d + 1.0;
  const mat& B0 = 2 * (kappa + 1) / kappa / c0 * Omega;
  // Rcout << c0 << endl << B0 << endl;
  
  // Main Loop
  for (int t = 1; t < T; t++) {
    //
    // Get new data
    vec yt1 = y.row(t).t();
    
    // Random number
    double u = unif_rand();
    
    // Probabilities of cluster
    vec dpredictive(mt + 1);
    dpredictive[0] =  alpha / (alpha + t) * dSt(y.row(1).t(), lambda, B0, c0);
    for (int i = 0; i < mt; i++) {
      vec atj = kappa * lambda + ntj[i] * mutj[i] / (kappa + ntj[i]);
      mat Dtj = ;
      double ctj = ;
      mat Btj = 2.0 * (kappa + ntj[i] + 1.0) / (kappa + ntj[0]) / ctj * (Omega  + 0.5 Dtj);
      dpredictive[i + 1] =  ntj[i] / (alpha + t) * dSt(y.row(1).t(), lambda, B0, c0);
    }


  }
  return T;
}

