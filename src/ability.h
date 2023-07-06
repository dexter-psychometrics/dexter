#ifndef DX_ITEMS_
#define DX_ITEMS_

#include <RcppArmadillo.h>

double Escore_single(double theta, const arma::vec& b, const arma::ivec& a, const arma::ivec& first,  const arma::ivec& last, const int n, const int max_a);

#endif