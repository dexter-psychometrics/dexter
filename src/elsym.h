#ifndef DEXTER_elsym_
#define DEXTER_elsym_

#include <RcppArmadillo.h>

template <class M, class V>
void elsym_i(const arma::vec& b_in, const arma::ivec& a_in, 
			int* ptr_first, int* ptr_last, 
			const int nit_all, const int exclude_item,
			M& g, V& g_all,
			const int from_item=0, const bool compute_full_elsym=true); 
					
					
template <class V>
void elsym(const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, V& g, V& gw, const int omit_item=-1);

#endif