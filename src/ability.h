#ifndef DEXTER_ABILITY_
#define DEXTER_ABILITY_

#include <RcppArmadillo.h>
/*
template<bool LOGLIK>
void deriv_theta(const double theta, const arma::vec& b, const arma::ivec& a, int* first, int* last, const int nit, const int max_a, arma::vec& PA, double& E, double& I, double &J);
*/

arma::vec ML_theta_c(const double score, const arma::mat& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last);
#endif