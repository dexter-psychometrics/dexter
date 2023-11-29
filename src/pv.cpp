#include <random>
#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include "myomp.h"
#include "priors.h"
#include "progress.h"

using namespace arma;

using Rcpp::Named;

#define SEED std::round(R::runif(0,1) * 2147483647)

// [[Rcpp::export]]
int omp_ncores()
{
	return omp_get_max_threads();
}

// to do: move to shared since it is useful elsewhere
ivec cum_iter_ivec(const ivec& v)
{
	const int n = v.n_elem;
	ivec out(n+1);
	out[0] = 0;
	out.tail(n) = cumsum(v);
	return out;
}


// [[Rcpp::export]]
Rcpp::List pv_chain_normal(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, const arma::ivec& first, const arma::ivec& last, 
			 const arma::ivec& bk_cnit, const arma::ivec& bk_max_a,
		     const arma::ivec& const_scoretab, const arma::ivec& scoretab_bk, const arma::ivec& scoretab_pop, 
			 const arma::ivec& scoretab_nscores, const arma::ivec& scoretab_np,
			 const arma::mat& mu_start, const arma::vec& sigma_start, const int npv, 
			 const arma::ivec progress_init, const int max_cores, const int warmup=10, const int step=1)
{
	const int ntab = scoretab_bk.n_elem;
	const int nchains = mu_start.n_cols;
	const int N = accu(scoretab_np);
	
	progress_prl pb(step*(npv-nchains) + nchains*warmup, progress_init);
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	
	ivec scoretab_cnp = cum_iter_ivec(scoretab_np);
	
	ivec scoretab_cnscores = cum_iter_ivec(scoretab_nscores);
	
	const ivec cscoretab = cumsum(const_scoretab); // not for iteration, so no zero added at the front	
	
	mat theta(N, std::max(nchains,npv)); 
	
	cube prior_log(mu_start.n_rows+3, (step*npv)/nchains + warmup + 1, nchains);
	
	int bstep = 0;
	if(bmat.n_cols > 1)
	{
		const int chain_iter = step * (npv/nchains)  + warmup;
		bstep = std::max(1u, bmat.n_cols/chain_iter);
	}
	

#pragma omp parallel num_threads(max_cores)
	{	
		const int thread = omp_get_thread_num();
		
		ivec scoretab(const_scoretab.n_elem); 
		vec mu(mu_start.n_rows);
		double sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
		
		int x,y,k, pvcol, cntr;
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
				scoretab = const_scoretab; // will get eaten
				
				mu = priors.theta;
				sigma = priors.sigma;	
				b = bmat.col(bcol);
				cntr = 0;
				
				for(int tab=0; tab < ntab; tab++)
				{
					if(pb.interrupted()) break;
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
						if(cntr++ > 50000)
						{
							if(thread==0) pb.checkInterrupt(); 
							if(pb.interrupted()) break;
							cntr = 0;
						}
					}	
				}
				
				prior_log.slice(chain).col(iter) = priors.as_vec();
				priors.update(lrng, theta.col(pvcol), scoretab_pop, scoretab_np, scoretab_cnp);
				bcol += bstep;
				if(bcol >= bmat.n_cols) bcol = chain;
				if(iter >= warmup && (iter - warmup) % step == 0) pvcol++;		
				pb.tick(thread == 0);	
			}
		}
	}
	pb.tick(true);
	if(pb.interrupted())
	{ 
		Rcpp::stop("user interrupt");
	}
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log);
}



// [[Rcpp::export]]
Rcpp::List pv_chain_mix(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, const arma::ivec& first, const arma::ivec& last, 
			 const arma::ivec& bk_cnit, const arma::ivec& bk_max_a,
		     const arma::ivec& gscoretab, const arma::ivec& gscoretab_bk, 
			 const arma::ivec& gscoretab_nscores, const arma::ivec& gscoretab_np,
			 const arma::mat& mu_start, const arma::mat& sigma_start, const arma::vec& p_start, const int npv, 
			 const arma::ivec progress_init, const int max_cores,
			 const int warmup=10, const int step=1)
{
	const int ntab = gscoretab_bk.n_elem;
	const int nchains = mu_start.n_cols;
	const int N = accu(gscoretab_np), len_gscoretab = gscoretab.n_elem;
	
	progress_prl pb(step*(npv-nchains) + nchains*warmup, progress_init);
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	
	mat theta(N, std::max(nchains,npv)); 
	
	cube prior_log(5, (step*npv)/nchains + warmup + 1, nchains);
	
	int bstep = 0;
	if(bmat.n_cols > 1)
	{
		const int chain_iter = step * (npv/nchains)  + warmup;
		bstep = std::max(1u, bmat.n_cols/chain_iter);	
	}
	
	ivec gscoretab_cnscores = cum_iter_ivec(gscoretab_nscores);
	
#pragma omp parallel num_threads(max_cores)
	{	
		const int thread = omp_get_thread_num();
		
		imat scoretab(gscoretab.n_elem, 2), cscoretab(gscoretab.n_elem, 2);		

		double mu, sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
	
		int x,y,k, pvcol, prs, cntr;
		double u, atheta;

		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);			
		
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
				// make two scoretab cols, one for each prior
				prs = 0;
				for(int i=0; i<len_gscoretab; i++)
					for(int sp=0; sp < gscoretab[i]; sp++)
						scoretab.at(i, priors.pop[prs++])++;
				
				// little trickery to get theta's in the right (= relating to booklet and score) position in output
				cscoretab.col(1) = cumsum(sum(scoretab,1));			
				cscoretab.col(0) = cscoretab.col(1) - scoretab.col(1);
				
				b = bmat.col(bcol);
				
				cntr = 0;
				for(int prior_num=0; prior_num <=1; prior_num++)
				{
					mu = priors.mu[prior_num];
					sigma = priors.sigma[prior_num];

					for(int tab=0; tab<ntab; tab++)
					{					
						if(pb.interrupted()) break;
						
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
							if(cntr++ > 50000)
							{
								if(thread==0) pb.checkInterrupt(); 
								if(pb.interrupted()) break;
								cntr = 0;
							}
						}						
					}
				}	
				prior_log.slice(chain).col(iter) = priors.as_vec();
				priors.upd_normal(lrng, theta.col(pvcol));
				bcol += bstep;
				if(bcol >= bmat.n_cols) bcol = chain;
				if(iter >= warmup && (iter - warmup) % step == 0) pvcol++;		
				pb.tick(thread == 0);	
			}
		}
	}
	
	pb.tick(true);
	
	if(pb.interrupted())
	{ 
		Rcpp::stop("user interrupt");
	}
	
	
	// The above will have sorted the theta's by booklet, booklet_score, assigned pop
	// the last part of the sort is undesirable since it may introduce an artefact or bias
	// so we shuffle within (booklet, booklet_score) for each column in the resulting pv's
	// tested that shuffle is correctly within bk-score and does actually shuffle
	
	// define gcscoretab
	ivec gcscoretab = cum_iter_ivec(gscoretab);	

	
// this seems to be seed safe, but might not be depending on what else is going on on the computer
#pragma omp parallel num_threads(max_cores)	
	{	
		const int thread = omp_get_thread_num();
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.jump(thread + 1);
		
#pragma omp for		
		for(int pvcol=0; pvcol < npv; pvcol++)
		{
			auto itr = theta.begin_col(pvcol);
			for(int i=0; i<len_gscoretab; i++)
			{
				if(gscoretab[i] > 1)
				{
					std::shuffle(itr + gcscoretab[i], itr + gcscoretab[i+1], lrng);
				}		
			}
		}
	}
	
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log);
}






// [[Rcpp::export]]
void PV_sve(const arma::vec& b, const arma::ivec& a, const arma::ivec& bk_first, const arma::ivec& bk_last, 					
			const arma::ivec& bcni,
			const arma::ivec& booklet_id, const arma::ivec& booklet_score, const arma::vec& mu, const double sigma, 
			const int max_cores,
			arma::mat& pv_mat, const int pv_col_indx=0, const int niter=1)
{
	const int np = pv_mat.n_rows;
	const int maxA = max(a);
	
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution prl_runif(0, 1);
	dqrng::normal_distribution prl_rnorm(0, 1);
	
	vec pv(pv_mat.colptr(pv_col_indx),np, false, true);
	
	
#pragma omp parallel num_threads(max_cores)	
	{
		const int thread = omp_get_thread_num();
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);
		
		double theta, u, acc;
		
		int x, bk, k;
		
		vec lookup(maxA+1);
		vec p(maxA+1, fill::zeros); 
		lookup[0] = 1.0;
	
#pragma omp for		
		for(int prs=0; prs<np; prs++)
		{
			bk = booklet_id[prs];
			for(int iter=0; iter<niter; iter++)
			{
				theta = prl_rnorm(lrng) * sigma + mu[prs];
				for(int j=1;j<=maxA;j++) 
					lookup[j] = std::exp(j*theta);
				
				x=0;
				for(int i=bcni[bk]; i<bcni[bk+1]; i++)
				{
					p[0] = b[bk_first[i]]; 
					k=1;
					for (int j=bk_first[i]+1;j<=bk_last[i];j++) 
					{
						p[k] = p[k-1] + b[j]*lookup[a[j]]; 
						k++;
					}
					u = p[k-1] * prl_runif(lrng);
					k = 0;
					while (u > p[k])
						k++;
					x += a[bk_first[i]+k];
				}
				
				acc = std::exp((theta-pv[prs])*(booklet_score[prs]-x));
				if(prl_runif(lrng)<acc)
					pv[prs] = theta;
			}
		}
	}
}

