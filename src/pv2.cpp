#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include <omp.h>
#include "priors.h"
#include "progress.h"

using namespace arma;

using Rcpp::Named;

#define SEED std::round(R::runif(0,1) * 2147483647)

// [[Rcpp::export]]
int ncores()
{
	return omp_get_max_threads();
}



// [[Rcpp::export]]
Rcpp::List pv_chain_normal(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, const arma::ivec& first, const arma::ivec& last, 
			 const arma::ivec& bk_cnit, const arma::ivec& bk_max_a,
		     const arma::ivec& const_scoretab, const arma::ivec& scoretab_bk, const arma::ivec& scoretab_pop, 
			 const arma::ivec& scoretab_nscores, const arma::ivec& scoretab_np,
			 const arma::mat& mu_start, const arma::vec& sigma_start, const int npv, 
			 const arma::ivec progress_init, const int warmup=10, const int step=1)
{
	const int ntab = scoretab_bk.n_elem;
	const int nchains = mu_start.n_cols;
	const int N = accu(scoretab_np);
	
	progress_prl pb(step*npv + nchains*warmup, progress_init);
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	

	
	ivec scoretab_cnp(ntab+1);
	scoretab_cnp[0] = 0;
	scoretab_cnp.tail(ntab) = cumsum(scoretab_np);
	
	ivec scoretab_cnscores(ntab+1);
	scoretab_cnscores[0] = 0;
	scoretab_cnscores.tail(ntab) = cumsum(scoretab_nscores);
	
	const ivec cscoretab = cumsum(const_scoretab); // not for iteration, so no zero added at the front	
	
	mat theta(N, std::max(nchains,npv)); 
	
	cube prior_log(mu_start.n_rows+3, (step*npv)/nchains + warmup + 1, nchains);
	
	int bstep = 0;
	if(bmat.n_cols > 1)
	{
		const int chain_iter = step * (npv/nchains)  + warmup;
		bstep = bmat.n_cols/chain_iter;	
	}
	

#pragma omp parallel
	{	
		const int thread = omp_get_thread_num();
		
		ivec scoretab(const_scoretab.n_elem); 
		vec mu(mu_start.n_rows);
		double sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
		
		int x,y,k, pvcol;
		double u, atheta;

		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);			
		
#pragma omp for
		for(int chain=0; chain<nchains; chain++)
		{
			hnorm_prior priors(mu_start.col(chain), sigma_start[chain]);
			
			pvcol = chain*(npv/nchains) + std::min(npv % nchains, chain);
			
			const int niter = (npv/nchains + (chain < npv % nchains) -1)*step + warmup;
			
			int bcol = 0;
			if(bstep>0) bcol = chain;
			
			for(int iter=0; iter<=niter; iter++)
			{	
				if(pb.interrupted) break;
				
				scoretab = const_scoretab; // will get eaten
				
				mu = priors.theta;
				sigma = priors.sigma;	
				b = bmat.col(bcol);

				for(int tab=0; tab < ntab; tab++)
				{
					int np = scoretab_np[tab];
					int max_score = scoretab_nscores[tab] - 1;			
					
					const int scoretab_start = scoretab_cnscores[tab];
					
					const int bk = scoretab_bk[tab];
					const int item_start = bk_cnit[bk];
					const int item_end = bk_cnit[bk+1];
					const int max_a = bk_max_a[bk];
					
					const double pmu = mu[scoretab_pop[tab]];
					
					expat[0] = 1.0;
					
					for(;max_score>=0; max_score--)
						if(scoretab[scoretab_start + max_score] > 0)
							break;
					
					while(np>0)
					{
						atheta = prl_rnorm(lrng) * sigma + pmu;
						
						for(int j=1; j<=max_a; j++) 
							expat[j] = std::exp(j*atheta);
						
						x = 0;
						for (int i=item_start; i<item_end; i++)
						{
							p[0] = b[first[i]]; 
							k=1;
							for (int j=first[i]+1;j<=last[i];j++) // note the +1
							{
								p[k] = p[k-1] + b[j] * expat[a[j]]; 
								k++;
							}
							u = p[k-1]*prl_runif(lrng);
							k=0;
							while (u > p[k])
								k++;
							if (k > 0) // if seems unnecessary since zero cat should be in a and A?
								x += A[first[i]+k];
								
							if(x > max_score)
								break;
						}	
						y = scoretab_start + x;
						
						if(scoretab[y] > 0)
						{
							theta.at(cscoretab[y] - scoretab[y], pvcol) = atheta;
							scoretab[y]--;
							np--;

							if(x == max_score && scoretab[y] == 0)
								for(;max_score>=0; max_score--)
									if(scoretab[scoretab_start + max_score] > 0)
										break;
						}		
					}	
				}
				
				prior_log.slice(chain).col(iter) = priors.as_vec();
				priors.update(lrng, theta.col(pvcol), scoretab_pop, scoretab_np, scoretab_cnp);
				bcol += bstep;							
				if(iter >= warmup && (iter - warmup) % step == 0) pvcol++;		
				pb.tick(thread == 0);	
			}
		}
	}
	if(pb.interrupted) Rcpp::stop("user interrupt");
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log);
}



// [[Rcpp::export]]
Rcpp::List pv_chain_mix(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, const arma::ivec& first, const arma::ivec& last, 
			 const arma::ivec& bk_cnit, const arma::ivec& bk_max_a,
		     const arma::ivec& gscoretab, const arma::ivec& gscoretab_bk, 
			 const arma::ivec& gscoretab_nscores, const arma::ivec& gscoretab_np,
			 const arma::mat& mu_start, const arma::mat& sigma_start, const arma::vec& p_start, const int npv, 
			 const arma::ivec progress_init,
			 const int warmup=10, const int step=1)
{
	const int ntab = gscoretab_bk.n_elem;
	const int nchains = mu_start.n_cols;
	const int N = accu(gscoretab_np), len_gscoretab = gscoretab.n_elem;
	
	progress_prl pb(step*npv + nchains*warmup, progress_init);
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	
	mat theta(N, std::max(nchains,npv)); 
	
	cube prior_log(5, (step*npv)/nchains + warmup + 1, nchains);
	
	int bstep = 0;
	if(bmat.n_cols > 1)
	{
		const int chain_iter = step * (npv/nchains)  + warmup;
		bstep = bmat.n_cols/chain_iter;	
	}
	
	ivec gscoretab_cnscores(ntab+1);
	gscoretab_cnscores[0] = 0;
	gscoretab_cnscores.tail(ntab) = cumsum(gscoretab_nscores);
	
#pragma omp parallel
	{	
		const int thread = omp_get_thread_num();
		
		imat scoretab(gscoretab.n_elem, 2), cscoretab(gscoretab.n_elem, 2);		

		double mu, sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
	
		int x,y,k, pvcol, prs;
		double u, atheta;

		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);			
		
		//int cntr=0;
#pragma omp for
		for(int chain=0; chain<nchains; chain++)
		{	
			mixture_prior priors(lrng, p_start(chain), mu_start.col(chain), sigma_start.col(chain), N);
			
			pvcol = chain*(npv/nchains) + std::min(npv % nchains, chain);
			
			const int niter = (npv/nchains + (chain < npv % nchains) -1)*step + warmup;
			
			int bcol = 0;
			if(bstep>0) bcol = chain;
			
			scoretab.zeros();
	
			for(int iter=0; iter<=niter; iter++)
			{
				if(pb.interrupted) break;
				// make two scoretab cols, one for each prior
				prs = 0;
				for(int i=0; i<len_gscoretab; i++)
					for(int sp=0; sp < gscoretab[i]; sp++)
						scoretab.at(i, priors.pop[prs++])++;
				
				// little trickery to get theta's in the right position in output
				// all pop 0 with same score will come before pop 1 theta's always though
				// this is not desirable
				// we could also randomize a bit
				// or we would need to keep track of old and new pop which is rather complicated
				// but that order matters is also not good
		
				cscoretab.col(1) = cumsum(sum(scoretab,1));			
				cscoretab.col(0) = cscoretab.col(1) - scoretab.col(1);
				
				b = bmat.col(bcol);
				
				for(int prior_num=0; prior_num <=1; prior_num++)
				{
					mu = priors.mu[prior_num];
					sigma = priors.sigma[prior_num];

					for(int tab=0; tab<ntab; tab++)
					{					
						const int scoretab_start = gscoretab_cnscores[tab]; // startpoint in scoretab 
						const int bk = gscoretab_bk[tab];
						
						const int item_start = bk_cnit[bk];
						const int item_end = bk_cnit[bk+1];
						const int max_a = bk_max_a[bk];
						
						expat[0] = 1.0;
						
						int max_score = gscoretab_nscores[tab] - 1;	
						
						for(;max_score>=0; max_score--)
							if(scoretab.at(scoretab_start + max_score, prior_num) > 0)
								break;
								
						int np = accu(scoretab.col(prior_num).subvec(scoretab_start, gscoretab_cnscores(tab+1)-1));				
						
						while(np>0)
						{
							atheta = prl_rnorm(lrng) * sigma + mu;
							
							for(int j=1; j<=max_a; j++) 
								expat[j] = std::exp(j*atheta);
							
							x = 0;
							for (int i=item_start; i<item_end; i++)
							{
								p[0] = b[first[i]]; 
								k=1;
								for (int j=first[i]+1; j<=last[i]; j++) // note the +1
								{
									p[k] = p[k-1] + b[j] * expat[a[j]]; 
									k++;
								}
								u = p[k-1] * prl_runif(lrng);
								k=0;
								while (u > p[k])
									k++;
								if (k > 0) 
									x += A[first[i]+k];
									
								if(x > max_score)
									break;
							}	
							y = scoretab_start + x;
							/*
							if(++cntr > 50000)
							{
								Rcpp::checkUserInterrupt();
								cntr = 0;
							}*/

							if(scoretab.at(y, prior_num) > 0)
							{
								theta.at(cscoretab.at(y, prior_num) - scoretab.at(y, prior_num), pvcol) = atheta;		
								
								scoretab.at(y, prior_num)--;
								np--;

								if(x == max_score && scoretab.at(y, prior_num) == 0)
									for(;max_score>=0; max_score--)
										if(scoretab.at(scoretab_start + max_score, prior_num) > 0)
											break;
							}	
						}						
					}
				}	
				prior_log.slice(chain).col(iter) = priors.as_vec();
				priors.upd_normal(lrng, theta.col(pvcol));
				bcol += bstep;
				if(iter >= warmup && (iter - warmup) % step == 0) pvcol++;		
				pb.tick(thread == 0);	
			}
		}
	}
	if(pb.interrupted) Rcpp::stop("user interrupt");
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log);
}

