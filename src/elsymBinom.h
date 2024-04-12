#ifndef DEXTER_lbnm_elsym_
#define DEXTER_lbnm_elsym_

#include <RcppArmadillo.h>

void elsym_i_binom(const arma::mat& lbinom,const arma::vec& b_in, const arma::ivec& a_in, 
			int* ptr_first, int* ptr_last, 
			const int nit_all, const int exclude_item,
			arma::mat& g, arma::vec& g_all,
			const int from_item=0, const bool compute_full_elsym=true); 


void elsym_binom(const arma::mat& lbinom, const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, arma::vec& g, arma::vec& gw, const int omit_item=-1);

#endif