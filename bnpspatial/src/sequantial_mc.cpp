#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES
#include <omp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
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
// samples an integer 0,...,n-1
inline int sample_int(vec prob) {
  vec probsum = cumsum(prob) / sum(prob);
  double u = unif_rand();
  int i;
  for (i = 0; u > probsum[i]; i++) {}
  return i;
}

// [[Rcpp::export]]
List dp_gaussian_mixture(
    const arma::mat& y,
    const double alpha,
    const arma::vec& lambda,
    const double kappa,
    const double nu,
    const arma::mat& Omega
    ) {
  
  // Dimensions of the problem
  uword T = y.n_rows;
  uword d = y.n_cols;
  
  // Initialization of essential state particles: no storage of past, T max number of clusters
  uword kt = 0; // Allocation at time t
  uword mt = 1; // Total components at time t,
  uvec ntj(T, fill::zeros);
  ntj[0] = 1; // Count per cluster
  mat meantj(d, T, fill::zeros); // Cluster mean
  meantj.col(0) = y.row(0).t();
  cube Stj(d, d, T, fill::zeros); // Cluster variance
  
  // Prior sampling constants
  double c0 = 2 * nu - d + 1.0;
  const mat& B0 = 2 * (kappa + 1) / kappa / c0 * Omega;
  
  // Main Loop
  for (uword t = 1; t < T; t++) {
    //
    // New data point
    vec ynew = y.row(t).t();
    
    // New cluster density
    vec dpredictive(mt + 1);
    
    for (uword i = 0; i < mt; i++) {
      // Cluster probabilities
      vec atj = (kappa * lambda + ntj[i] * meantj.col(i)) / (kappa + ntj[i]);
      mat Dtj = Stj.slice(i) + kappa * ntj[i] / 
        (kappa + ntj[i]) * (lambda - meantj.col(i)) * (lambda - meantj.col(i)).t();
      double ctj = 2 * nu + ntj[i] - d + 1.0;
      mat Btj = 2.0 * (kappa + ntj[i] + 1.0) / (kappa + ntj[i]) / ctj * (Omega  + 0.5 * Dtj);
      // Rcout << "ynew" << ynew << "\nntj" << ntj[i] << "\nBtj" << Btj << "\nDtj" << Dtj << "\natj " << atj << "\nctj" << ctj << endl;
      dpredictive[i] =  ntj[i] / (alpha + t) * dSt(ynew, atj, Btj, ctj);
    }
    
    // Prior probability
    dpredictive[mt] =  alpha / (alpha + t) * dSt(ynew, lambda, B0, c0);
    
    // Update values
    kt = sample_int(dpredictive);
    // Rcout << "dpredictive\n" << dpredictive << "kt: " << kt << "\nmt: " << mt << "\nntjkt: " << ntj[kt] << endl;
    kt == mt && mt++;
    ntj[kt] += 1;
    vec meantjt = meantj.col(kt);
    meantj.col(kt) =  ((ntj[kt] - 1) * meantj.col(kt) + ynew) / ntj[kt];
    Stj.slice(kt) += ynew * ynew.t() + 
      (ntj[kt] - 1) * meantjt * meantjt.t() - ntj[kt] * meantj.col(kt) * meantj.col(kt).t();
  }
  
  // Output List
  ntj.resize(mt);
  meantj.resize(d, mt);
  Stj.resize(d, d, mt);
  return Rcpp::List::create(
    Named("m") = mt,
    Named("nj") = ntj,
    Named("meanj") = meantj,
    Named("Sj") = Stj);
}

