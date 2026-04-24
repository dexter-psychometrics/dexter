#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include "myomp.h"
#include "shared.h"
#include "ability.h"

#define SEED std::round(R::runif(0,1) * 2147483647)

using namespace arma;

/*
imputation of plausible scores

b has either one column or nbr of columns equal to pv

first, last denote the plausible item set

person_id, item_first, item_score are columns denoting existing scores
ordered by person_id, item_first
item_first == first means the item was already answered and score should be kept
item_first contains only items that exist in first
person_id contains 0 based indexes of pv, terminated by -1 (person_id can therefore be one element longer than item_id and item_score) 
*/
template <bool by_item>
arma::imat impute_NRM_tpl(const arma::mat& pv, const arma::mat& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last,
							 const arma::ivec& person_id, const arma::ivec& item_first, const arma::ivec& item_score,
							const int max_cores)
{
	const bool multiple_b = b.n_cols>1;
	const int npv = pv.n_cols, np = pv.n_rows, nit = first.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	
	int out_n_rows = pv.n_rows;
	if(by_item) out_n_rows *= nit;
	imat score(out_n_rows, pv.n_cols, fill::zeros);		
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution runif(0, 1);
	
#pragma omp parallel num_threads(max_cores)
{
	const int thread = omp_get_thread_num();
	vec p(maxA+2, fill::zeros);
	p[0] = 1;
	vec expat(maxA+1);
	expat[0] = 1;
	
	bool any_obs = false;
	int out_indx, k, indx_obs, bcol=0;
	double u;
	
	dqrng::xoshiro256plus lrng(rng);      		
	lrng.long_jump(thread + 1);	
	
#pragma omp for	
	for(int pvcol=0; pvcol < npv; pvcol++)
	{
		if(multiple_b) bcol = pvcol;
		out_indx = 0;
		indx_obs = 0;
		for(int prs=0; prs < np; prs++)
		{
			any_obs = (prs == person_id[indx_obs]);

			for(int i=1;i<=maxA;i++){ expat[i] = std::exp(i*pv.at(prs, pvcol));}
			
			for (int i=0; i<nit; i++)
			{
				if(any_obs && item_first[indx_obs] == first[i])
				{
					if(by_item) score.at(out_indx++, pvcol) = item_score[indx_obs++];
					else score.at(prs, pvcol) += item_score[indx_obs++];
					any_obs = (prs == person_id[indx_obs]);
					continue;
				}
				
				k=1;
				for (int j=first[i];j<=last[i];j++) 
				{
					p[k] = p[k-1] + b.at(j, bcol)*expat[a[j]];
					k++;
				}
				u = p[k-1] * runif(lrng);
				k = 0;
				while (u>p[k]) {k++;}
				if(k>0)
				{
					if (by_item) score.at(out_indx, pvcol) = a[first[i]+k-1];
					else score.at(prs, pvcol) += a[first[i]+k-1];
				}
				if(by_item) out_indx++;
			}
		}
	}
}
	return score;
}


//[[Rcpp::export]]
arma::imat impute_NRM_C(const arma::mat& pv, const arma::mat& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last,
							 const arma::ivec& person_id, const arma::ivec& item_first, const arma::ivec& item_score,
							const bool by_item, const int max_cores)
{
	if(by_item)
	{
		return impute_NRM_tpl<true>(pv,b, a, first, last,	
										person_id, item_first, item_score, max_cores);
	} else
	{
		return impute_NRM_tpl<false>(pv,b, a, first, last,	
										person_id, item_first, item_score, max_cores);
	}
}


// combining this with previous will not work, since the existing value cannot be combined with parallel by person, only by pvcol
// we could combine them in the cpp wrapper (but just as easy in R wrapper)

template <bool by_item>
arma::imat sample_NRM_tpl(const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
						  const arma::ivec& first, const arma::ivec& last, const int max_cores)
{

	const int nit = first.n_elem, np = theta.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int ncol = by_item ? nit : 1;
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution runif(0, 1);
  
	imat score(np, ncol, fill::zeros);  
  
#pragma omp parallel num_threads(max_cores)
{
	const int thread = omp_get_thread_num();
	int k;
	double u;
	vec p(maxA+3);
	vec expat(maxA+1);
	expat[0] = 1;
	p[0] = 1;
	
	dqrng::xoshiro256plus lrng(rng);      		
	lrng.long_jump(thread + 1);	
	
#pragma omp for	
	for (int prs=0; prs < np; prs++)
	{
		for(int i=1; i<=maxA; i++)
			expat[i] = std::exp(i*theta[prs]);
    
		for (int i=0; i<nit; i++)
		{
			k = 1;
			for (int j=first[i]; j<=last[i]; j++)
			{
				p[k] = p[k-1]+b[j]*expat[a[j]];
				k++;
			}
			u = p[k-1] * runif(lrng);
			k = 0;
			while(u>p[k]) k++;
			
			if(by_item && k > 0) score.at(prs, i) = a[first[i]+k-1];
			if(!by_item && k > 0) score[prs] += a[first[i]+k-1];			
		}
	}
}
	return score;
}

//[[Rcpp::export]]
arma::imat sample_NRM_C(const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
						  const arma::ivec& first, const arma::ivec& last, const bool by_item, const int max_cores)
{
	if(by_item) return sample_NRM_tpl<true>(theta, b, a, first, last, max_cores);
	else return sample_NRM_tpl<false>(theta, b, a, first, last, max_cores);
}



// [[Rcpp::export]]
arma::imat sampleIMC(const arma::vec& bIM, const arma::vec& cIM, const arma::ivec& a, arma::ivec& first, arma::ivec& last,
					const arma::ivec& scoretab)
{
	// scoretab must already be based on the sample and must include 0 scores even for impossible
	const int nit = first.n_elem;
	const int max_score = scoretab.n_elem - 1;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int nP = accu(scoretab);
	
	double u, theta = -2;

	int score, k, pi;
	
	
	
	vec wmem(maxA+1);
	vec logb = log(bIM), logc = log(cIM);
	
	vec b(&bIM[0], bIM.n_elem);

	ivec cs_scoretab(max_score+1);
	vec lookup(maxA+1);
	vec p(maxA+3);
	
	imat out(nP, nit);	
	
	lookup[0] = 1;
	p[0] = 1;
	
	cs_scoretab[0] = 0;
	for(int i=0; i<max_score; i++)
		cs_scoretab[i+1] = scoretab[i] + cs_scoretab[i];

	
	for(int s=1; s<max_score; s++)
	{
		for(int i=0; i<nit; i++)
			for (int j=first[i]; j<=last[i]; j++) 
				b[j] = exp(logb[j] + s * a[j] * logc[i]);
		
		theta = as_scalar(ML_theta_c((double)s, b, a, first, last));
		
		for(int i=1; i<=maxA; i++)
			lookup[i] = std::exp(i*theta);
				
		pi = cs_scoretab[s];
		// sample
		while(pi < cs_scoretab[s+1])
		{
			score = 0;
			out.row(pi).zeros();
			for (int i=0;i<nit;i++)
			{				
				k = 1;				
				for (int j=first[i]; j<=last[i]; j++) 
				{
					p[k] = p[k-1]+b[j]*lookup[a[j]];
					k++;
				}
				u = p[k-1]*R::runif(0,1);
				k = 0;
				while(u>p[k]) k++;
				if(k>0) 
				{
					score += a[first[i]+k-1];
					out.at(pi, i) = a[first[i]+k-1];
					
					if(score > s ) break;
				}				
			}
			if(score == s)
				pi++;	
		}
		theta += 0.1;
	}

	return out;
}

