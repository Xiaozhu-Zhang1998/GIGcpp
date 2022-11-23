// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// Algorithm 0: Sample from truncated exponential distribution
Eigen::VectorXd Rcpp_trunc_exp(int n, double rate, double a, double b) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  for (auto i = 0; i < n; ++i)
    x[i] = R::rexp(1/rate) + a;
  if (isinf(b))
    return x;
  double c = b - a;
  for (auto i = 0; i < n; ++i)
    x[i] = x[i] - floor((x[i] - a) / c) * c;
  return x;
}


// Algorithm 1: Sample from truncated gamma distribution
Eigen::VectorXd Rcpp_trunc_gamma(int n, double shape, double rate, double lwr) {
  double p = R::pgamma(lwr, shape, 1/rate, 0, 1);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
  for (auto i = 0; i < n; ++i)
    v[i] = R::qgamma(-R::rexp(1) + p, shape, 1/rate, 0, 1);
  return v;
}


// Algorithm 2: Rejection sampling
Eigen::VectorXd Rcpp_rejection_sampling(int N, double v, double b, Eigen::VectorXd cutpoint) {
  int len_cutpoint = cutpoint.size();
  Eigen::VectorXd ext_cutpoint = Eigen::VectorXd::Zero(len_cutpoint + 2);
  for (auto i = 0; i < len_cutpoint; ++i) {
    ext_cutpoint[i + 1] = cutpoint[i];
  }
  ext_cutpoint[len_cutpoint + 1] = 1.0 / 0;
  Eigen::VectorXd M = Eigen::VectorXd::Zero(len_cutpoint + 1);
  Eigen::VectorXd prob = Eigen::VectorXd::Zero(len_cutpoint + 1);
  double prob_sum = 0;
  for (auto i = 0; i < len_cutpoint; ++i) {
    M[i] = R::pgamma(1/cutpoint[i], -v, 2/b, 0, 0);
    if (i > 0) {
      prob[i] = R::pexp(cutpoint[i], 2/b, 1, 0) - R::pexp(cutpoint[i - 1], 2/b, 1, 0);
      prob[i] *= M[i];
      prob_sum += prob[i];
    }
  }
  M[len_cutpoint] = 1;
  prob[0] =  R::pexp(cutpoint[0], 2/b, 1, 0) * M[0];
  prob[len_cutpoint] = R::pexp(cutpoint[len_cutpoint - 1], 2/b, 0, 0);
  prob_sum += prob[0] + prob[len_cutpoint];
  
  Eigen::VectorXd Y = Eigen::VectorXd::Zero(N);
  int len = 0, cnt = 0;
  while (len < N) {
    double unif = R::runif(0, prob_sum);
    double tmp = unif;
    int ind = 0;
    while(ind <= len_cutpoint) {
      if (unif < prob[ind]) break;
      unif -= prob[ind++];
    }
    if (ind > len_cutpoint) ind = len_cutpoint;
    double y = Rcpp_trunc_exp(1, b/2, ext_cutpoint[ind], ext_cutpoint[ind + 1])[0];
    double rt = R::pgamma(1/y, -v, 2/b, 0, 0) / M[ind];
    ++cnt;
    if(R::runif(0, 1) <= rt) {
      Y[len++] = y;
    }
  }
  return Y;
}


// Algorithm 3: Find cutoff points given rejection rate
Eigen::VectorXd Rcpp_find_cutoff_under_rej_rate(double v, double b, double eps) {
  Eigen::VectorXd cutoff_points = Eigen::VectorXd::Zero(10);
  int cur_size = 0;
  double cur_threshold = 1 - eps;
  double left_area = 1, right_area = 0;
  double last_cutoff = 1.0 / 0;
  while (left_area > (left_area + right_area) * eps) {
    double new_cutoff = 1 / R::qgamma(cur_threshold, -v, 2/b, 0, 0);
    cutoff_points[cur_size++] = new_cutoff;
    if (cur_size == cutoff_points.size()) {
      Eigen::VectorXd new_cutoff_points = Eigen::VectorXd::Zero(2 * cutoff_points.size());
      for (int i = 0; i < cur_size; ++i)
        new_cutoff_points[i] = cutoff_points[i];
      cutoff_points = new_cutoff_points;
    }
    cur_threshold *= 1 - eps;
    double left_prop = R::pexp(new_cutoff, 2/b, 1, 0) / R::pexp(last_cutoff, 2/b, 1, 0);
    right_area += (1 - left_prop) * left_area;
    left_area *= left_prop * (1 - eps);
    last_cutoff = new_cutoff;
  }
  return cutoff_points.segment(0, cur_size).reverse();
}


// Algorithm 4: Find cutoff points given the desired number
Eigen::VectorXd Rcpp_find_cutoff_under_fixed_number(double v, double b, int n) {
  double l = 0, r = 1;
  while (r - l > 1e-6) {
    double m = (r + l) / 2;
    int n_m = Rcpp_find_cutoff_under_rej_rate(v, b, m).size();
    if (n_m < n) 
      r = m;
    else
      l = m;
  }
  return (Rcpp_find_cutoff_under_rej_rate(v, b, l));
}


// Algorithm 5: Sample from GIG

//' @title A Generator for Generalized Gaussian Distribution
//' @description A function that generates random variates from generalized Gaussian distributions. 
//' User can specify a rejection rate upper bound or a desired number of cutoff points.
//' @param N number of variates.
//' @param lambda parameter lambda.
//' @param psi parameter psi.
//' @param chi parameter chi.
//' @param eps desired rejection rate between 0 and 1; active only when `K` == 0. 
//' @param K desired number of cutoff points
//' @note See Zhang and Reiter (2022)
//' @return a vector of `N` random variates of the specified GIG distribution
// [[Rcpp::export]]

Eigen::VectorXd rgig(int N, double lambda, double psi, double chi, double eps = 0.5, int K = 0) {
  if(lambda == 0) {
    Rcpp::stop("lambda cannot be 0.");
  }
  if(eps <= 0 || eps >= 1) {
    Rcpp::stop("eps must be between 0 and 1.");
  }
  if(K < 0) {
    Rcpp::stop("K must be non-negative.");
  }
    
  bool flag = false;
  Eigen::VectorXd cutpoint;
  if(lambda > 0) {
    lambda = -lambda;
    std::swap(psi, chi);
    flag = true;
  }
  double alpha = std::sqrt(psi / chi);
  double beta = std::sqrt(psi * chi);
  if(K == 0)
    cutpoint = Rcpp_find_cutoff_under_rej_rate(lambda, beta, eps);
  else
    cutpoint = Rcpp_find_cutoff_under_fixed_number(lambda, beta, K);
  
  Eigen::VectorXd Y = Rcpp_rejection_sampling(N, lambda, beta, cutpoint);
  Eigen::VectorXd X = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < N; ++i) {
    X[i] = Rcpp_trunc_gamma(1, -lambda, beta/2, 1/Y[i])[0];
    X[i] = 1 / (alpha * X[i]);
    if(flag)
      X[i] = 1 / X[i];
  }
  return (X);
}


