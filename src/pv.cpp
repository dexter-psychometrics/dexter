#include <random>
#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include "myomp.h"
#include "elsym.h"
#include "priors.h"
#include "progress.h"
#include "shared.h"

using namespace arma;

using Rcpp::Named;

#define SEED std::round(R::runif(0,1) * 2147483647)

// [[Rcpp::export]]
int omp_ncores()
{
	return omp_get_max_threads();
}


// probabilities of sumscores based on log elsym, used for metropolis step

void px_theta(vec& px, const arma::vec& lg, const double theta, const bool log_prob = true)
{
	for(int score=0; score < px.n_elem; score++)
		px[score] = lg[score] + score * theta;
	
	px -= logsumexp(px);
	
	if(!log_prob)
		px = exp(px);
}

// using pointer, return only for single x
// not psbl to make any faster
double log_px_theta(const double *lg, const int max_score, const int x, const double theta)
{
	double m = 0, s;
	
	for(int score=1; score <= max_score; score++)
		m = std::max(*(lg+score) + score * theta, m);		

	s = std::exp(-m);
	for(int score=1; score <= max_score; score++)
		s += std::exp(*(lg+score) + score * theta-m);
	
	return *(lg+x) + x * theta - m - std::log(s);
}

/********** Plausible value algorithms with individual priors ************/


/*
Single value exchange
 - fast
 - works well for shorter booklets
 - high autocorrelation/low coverage for scores that are unlikley given the prior
 - adjusted to work with incomplete data (inserts draw from the prior)
 
*/

// [[Rcpp::export]]
void PV_sve(const arma::vec& b, const arma::ivec& a, const arma::ivec& bk_first, const arma::ivec& bk_last, 					
			const arma::ivec& bcni,
			const arma::ivec& booklet_id, const arma::ivec& booklet_score, const arma::vec& mu, const double sigma, 
			const int max_cores,
			arma::mat& pv_mat, const arma::ivec& missing_data, const int pv_col_indx=0, const int niter=1)
{
	const int np = pv_mat.n_rows;
	const int maxA = max(a);
	const bool complete = all(missing_data==0);
	
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution runif(0, 1);
	dqrng::normal_distribution rnorm(0, 1);
	
	vec pv(pv_mat.colptr(pv_col_indx),np, false, true);
	
	const ivec data_offset = cumsum(missing_data);
	
#pragma omp parallel num_threads(max_cores)	
	{
		const int thread = omp_get_thread_num();
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);
		
		double theta, u, acc;
		
		int x, bk, k, bk_score;
		
		vec exp_at(maxA+1);
		vec p(maxA+1, fill::zeros); 
		exp_at[0] = 1;
		p[0]=1;
#pragma omp for		
		for(int prs=0; prs<np; prs++)
		{
			if(complete)
			{
				bk = booklet_id[prs];
				bk_score = booklet_score[prs];
			}
			else
			{
				if(missing_data[prs] == 1)
				{				
					pv[prs] = rnorm(lrng) * sigma + mu[prs];
					continue;
				}	
				bk = booklet_id[prs-data_offset[prs]];
				bk_score = booklet_score[prs-data_offset[prs]];
			}
			
			for(int iter=0; iter<niter; iter++)
			{
				theta = rnorm(lrng) * sigma + mu[prs];
				for(int j=1; j<=maxA; j++) 
					exp_at[j] = std::exp(j*theta);
				
				x=0;
				for(int i=bcni[bk]; i<bcni[bk+1]; i++)
				{
					k=1;
					for (int j=bk_first[i];j<=bk_last[i];j++) 
					{
						p[k] = p[k-1] + b[j]*exp_at[a[j]]; 
						k++;
					}
					u = p[k-1] * runif(lrng);
					k = 0;
					while (u > p[k])
						k++;
					if(k>0) x += a[bk_first[i]+k-1];
				}
				acc = std::exp((theta-pv[prs])*(bk_score-x));				
				
				if(runif(lrng) < acc)
					pv[prs] = theta;
			}
		}
	}
}



/* 
Plausible Values with individual priors

Uses a simple metropolis step. Relatively slow but sure algorithm. Works well for cases that are pathological in the 2 other algorithms, i.e. unlikely scores in long booklets.

-  booklet_id, booklet_score and mu are by person
-  pv_res contains starting values/current pv's which will be overwritten

*/


// [[Rcpp::export]]
void pv_metro(const arma::vec& lg, const arma::ivec& booklet_maxscore, 
					const arma::ivec& booklet_id, const arma::ivec& booklet_score, const arma::vec& prior_mu, const double prior_sigma,
					arma::vec& pv_res, const int max_cores,  const int n_updates=10)
{
	const int np = booklet_score.n_elem;
	const double jump_sd = .8 * prior_sigma;	
	
	const ivec bk_cnscores = ivec2_iter(booklet_maxscore + 1);
	
	const int max_bkscore = max(booklet_maxscore);
	
	dqrng::xoshiro256plus rng(SEED);
	dqrng::uniform_distribution runif(0, 1);
	dqrng::normal_distribution rnorm(0, 1);
	
#pragma omp parallel num_threads(max_cores)
{
	const int thread = omp_get_thread_num();
	int bk;
	double pv, pv_new, p, p_new, acc;
	
	dqrng::xoshiro256plus lrng(rng);      		
	lrng.long_jump(thread + 1);	
	
#pragma omp for	
	for(int prs = 0; prs < np; prs++)
	{
		bk = booklet_id[prs];
		pv = pv_res[prs];
		p = log_px_theta(lg.memptr() + bk_cnscores[bk], booklet_maxscore[bk], booklet_score[prs], pv) + log_normpdf(pv, prior_mu[prs], prior_sigma); 
			
		for(int iter=0; iter < n_updates; iter++)
		{
			pv_new = rnorm(lrng) * jump_sd + pv; 
			p_new = log_px_theta(lg.memptr() + bk_cnscores[bk], booklet_maxscore[bk], booklet_score[prs], pv_new) + log_normpdf(pv_new, prior_mu[prs], prior_sigma); 
			acc = std::exp(p_new - p);
			if(runif(lrng) < acc)
			{
				pv = pv_new;
				p = p_new;
			}
		}			
		pv_res[prs] = pv;
	}
}
}


/********** Plausible value algorithms with group priors ************/


/* 
Same algorithm as for individual pv's, now using a scoretab and a common prior
-  starting value is a draw from the eap
*/

arma::vec pvmetro_stab(int* scoretab, const arma::vec& lg, const double prior_mu, const double prior_sigma, dqrng::xoshiro256plus& lrng, const int n_updates=50)
{
	const int max_score = lg.n_elem-1, ngrid=100;
	
	int np = 0;
	for(int score=0; score <= max_score; score++) np += *(scoretab + score);
	
	double post_mu, post_sigma, pv, pv_new, p, p_new, acc;
	
	dqrng::uniform_distribution runif(0, 1);
	dqrng::normal_distribution rnorm(0, 1);
	
	vec px(max_score+1), out(np), post(ngrid);
	
	vec grid = linspace(-6, 6, ngrid);
	vec w = normpdf(grid);
	w /= sum(w);
	grid = grid * prior_sigma + prior_mu;
	
	mat p_score(max_score+1, ngrid);
	for(int i=0; i < ngrid; i++)
	{
		px_theta(px, lg, grid[i], false);		
		p_score.col(i) = px;			
	}

	int prs = 0;
	for(int score=0; score <= max_score; score++) if(*(scoretab + score) > 0) 
	{
		// approximate mean and sigma of posterior P(theta|x+) = L*prior
		for(int i=0; i < ngrid; i++)
			post[i] = w[i] * p_score.at(score,i);

		post /= accu(post);
		post_mu = accu(grid % post);
		post_sigma = std::sqrt(accu(post % square(grid - post_mu)));

		// using only log probabilities
		for(int j=0; j<*(scoretab + score); j++)
		{
			pv = rnorm(lrng) * post_sigma + post_mu; // start
			px_theta(px, lg, pv);
			p = px[score] + log_normpdf(pv, prior_mu, prior_sigma); 
			for(int iter=0; iter < n_updates; iter++)
			{
				pv_new = rnorm(lrng) * post_sigma + pv; // jumping kernel
				px_theta(px, lg, pv_new);
				p_new = px[score] + log_normpdf(pv_new, prior_mu, prior_sigma); 
				acc = std::exp(p_new - p);
				if(runif(lrng) < acc)
				{
					pv = pv_new;
					p = p_new;
				}
			}
			out[prs++] = pv;		
		}
	}
	return out;
}


/*
Fill remaining scoretab when rejection algorithm becomes too inefficient
-  computes and overwrites log_gamma if not previously computed
*/

arma::vec draw_unlikely_pvs(arma::vec& lg, int *first, int *last, const int nit, const vec& b, const ivec& a, 
					int* scoretab, const int n_scores, const double prior_mu, const double prior_sigma, dqrng::xoshiro256plus& lrng)
{
	if(lg.n_elem == 0)
	{
		//lg = vec(n_scores);
		lg.set_size(n_scores);
		vec gw(n_scores);
		elsym(b, a, first, last, nit, lg, gw, -1);
		lg = log(lg);
	}
	return pvmetro_stab(scoretab, lg, prior_mu, prior_sigma, lrng);
}



/* 
Pv's using Rejection algorithm with 
   1) (hierarchical) normal prior
   2) mixture prior 
-  based on scoretabs   
-  pathological cases are resolved with draw_unlikely_pvs

-> to do: unlikley pv's can only use a, while the normal algorithm can use A 
	(not urgent since currently the public api guarantees that A is always equal to a)
*/
// [[Rcpp::export]]
Rcpp::List pv_chain_normal(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, arma::ivec& first, arma::ivec& last, 
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
	
	ivec scoretab_cnp = ivec2_iter(scoretab_np);
	
	ivec scoretab_cnscores = ivec2_iter(scoretab_nscores);
	
	const ivec cscoretab = cumsum(const_scoretab); // not for iteration, so no zero added at the front	
	
	mat theta(N, std::max(nchains,npv)); 
	
	cube prior_log(mu_start.n_rows+3, (step*npv)/nchains + warmup + 1, nchains);
	
	int bstep = 0;
	if(bmat.n_cols > 1)
	{
		const int chain_iter = step * (npv/nchains)  + warmup;
		bstep = std::max(1u, bmat.n_cols/chain_iter);
	}
	
	int n_alt=0;
#pragma omp parallel num_threads(max_cores)
	{	
		const int thread = omp_get_thread_num();
		
		ivec scoretab(const_scoretab.n_elem); 
		vec mu(mu_start.n_rows);
		double sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
		vec pv_add;
		
		int x,y,k, pvcol, cntr, np, np_prev;
		double u, atheta;
		
		expat[0] = 1;
		p[0]=1;
		
		field<vec> log_gamma(bk_max_a.n_elem);

		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);			
		
#pragma omp for reduction(+: n_alt)
		for(int chain=0; chain<nchains; chain++)
		{
			hnorm_prior priors(mu_start.col(chain), sigma_start[chain]);
			
			pvcol = chain*(npv/nchains) + std::min(npv % nchains, chain);
			
			const int niter = (npv/nchains + (chain < npv % nchains) -1)*step + warmup;
			
			int bcol = 0;
			if(bstep > 0) bcol = chain;

				
			for(int iter=0; iter<=niter; iter++)
			{					
				if(thread==0) pb.checkInterrupt(); 
				if(pb.interrupted()) break;
				
				scoretab = const_scoretab; // will get eaten
				
				mu = priors.theta;
				sigma = priors.sigma;	
				b = bmat.col(bcol);
				if(bstep>0) for(int bk=0; bk < log_gamma.n_elem; bk++) log_gamma(bk).reset(); // clear cache of log gamma
				cntr = 0;
				
				for(int tab=0; tab < ntab; tab++)
				{					
					np = scoretab_np[tab];
					np_prev = np;
					
					const int scoretab_start = scoretab_cnscores[tab];
					
					const int bk = scoretab_bk[tab];
					const int item_start = bk_cnit[bk];
					const int item_end = bk_cnit[bk+1];
					const int max_a = bk_max_a[bk];
					
					const double pmu = mu[scoretab_pop[tab]];				
					
					while(np>0)
					{
						atheta = prl_rnorm(lrng) * sigma + pmu;
						
						for(int j=1; j<=max_a; j++) 
							expat[j] = std::exp(j*atheta);
						
						x = 0;
						for (int i=item_start; i<item_end; i++)
						{
							k=1;
							for (int j=first[i];j<=last[i];j++) 
							{
								p[k] = p[k-1] + b[j] * expat[a[j]]; 
								k++;
							}
							u = p[k-1]*prl_runif(lrng);
							k=0;
							while (u > p[k])
								k++;
							if (k > 0) 
							{
								x += A[first[i]+k-1];								
							}	
							
						}	
						y = scoretab_start + x;
						
						if(scoretab[y] > 0)
						{
							theta.at(cscoretab[y] - scoretab[y], pvcol) = atheta;
							scoretab[y]--;
							np--;
						}
						if(cntr++ > 5000) // fail rejection sampling at >5000 consecutive draws without any new pv
						{
							if(np_prev == np)
							{
								pv_add = draw_unlikely_pvs(log_gamma(bk), first.memptr() + item_start, last.memptr() + item_start, item_end - item_start, b, a, 
																scoretab.memptr() + scoretab_start, scoretab_nscores[bk], pmu, sigma, lrng);
								n_alt += pv_add.n_elem;
								y = scoretab_start;
								for(int i=0; i < pv_add.n_elem; i++)
								{
									while(scoretab[y]==0) y++;
									theta.at(cscoretab[y] - scoretab[y], pvcol) = pv_add[i];
									scoretab[y]--;
								}
								break;
							}
							cntr = 0;
							np_prev = np;
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
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log, Named("n_alt_pv") = n_alt);
}

// mixture prior
// [[Rcpp::export]]
Rcpp::List pv_chain_mix(const arma::mat& bmat, const arma::ivec& a, const arma::ivec& A, arma::ivec& first, arma::ivec& last, 
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
	
	ivec gscoretab_cnscores = ivec2_iter(gscoretab_nscores);
	
#pragma omp parallel num_threads(max_cores)
	{	
		const int thread = omp_get_thread_num();
		
		imat scoretab(gscoretab.n_elem, 2), cscoretab(gscoretab.n_elem, 2);		

		double mu, sigma;
			
		vec b(bmat.n_rows), expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros), pv_add;
		field<vec> log_gamma(bk_max_a.n_elem);
	
		int x,y,k, pvcol, prs, cntr, np, np_prev, bcol;
		double u, atheta;
		
		expat[0] = 1;
		p[0] = 1;

		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);			
		
#pragma omp for
		for(int chain=0; chain<nchains; chain++)
		{	
			mixture_prior priors(lrng, p_start(chain), mu_start.col(chain), sigma_start.col(chain), N);
			
			pvcol = chain*(npv/nchains) + std::min(npv % nchains, chain);
			
			const int niter = (npv/nchains + (chain < npv % nchains) -1)*step + warmup;
			
			bcol = 0;
			if(bstep>0) bcol = chain;
			
			scoretab.zeros();			
	
			for(int iter=0; iter<=niter; iter++)
			{				
				if(pb.interrupted()) break;
				// make two scoretab cols, one for each prior
				prs = 0;
				for(int i=0; i<len_gscoretab; i++)
					for(int sp=0; sp < gscoretab[i]; sp++)
						scoretab.at(i, priors.pop[prs++])++;
				
				// little trickery to get theta's in the right (= relating to booklet and score) position in output
				cscoretab.col(1) = cumsum(sum(scoretab,1));			
				cscoretab.col(0) = cscoretab.col(1) - scoretab.col(1);
				
				b = bmat.col(bcol);
				if(bstep > 0) for(int bk=0; bk < log_gamma.n_elem; bk++) log_gamma(bk).reset(); // clear cache of log gamma
				cntr = 0;
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
						
								
						np = accu(scoretab.col(prior_num).subvec(scoretab_start, gscoretab_cnscores(tab+1)-1));				
						np_prev = np;
						while(np>0)
						{
							atheta = prl_rnorm(lrng) * sigma + mu;
							
							for(int j=1; j<=max_a; j++) 
								expat[j] = std::exp(j*atheta);
							
							x = 0;
							for (int i=item_start; i<item_end; i++)
							{
								k=1;
								for (int j=first[i]; j<=last[i]; j++) 
								{
									p[k] = p[k-1] + b[j] * expat[a[j]]; 
									k++;
								}
								u = p[k-1] * prl_runif(lrng);
								k=0;
								while (u > p[k])
									k++;
								if (k > 0) 
									x += A[first[i]+k-1];
							}	
							y = scoretab_start + x;

							if(scoretab.at(y, prior_num) > 0)
							{
								theta.at(cscoretab.at(y, prior_num) - scoretab.at(y, prior_num), pvcol) = atheta;		
								
								scoretab.at(y, prior_num)--;
								np--;
							}
							if(cntr++ > 5000)
							{
								if(np_prev == np)
								{
									pv_add = draw_unlikely_pvs(log_gamma(bk), first.memptr() + item_start, last.memptr() + item_start, item_end - item_start, b, a, 
																scoretab.colptr(prior_num) + scoretab_start, gscoretab_nscores[tab], mu, sigma, lrng);
	
									y = scoretab_start;
									for(int i=0; i < pv_add.n_elem; i++)
									{
										while(scoretab.at(y, prior_num)==0) y++;
										theta.at(cscoretab.at(y, prior_num) - scoretab.at(y, prior_num), pvcol) = pv_add[i];	
										scoretab.at(y, prior_num)--;
									}
									break;
								}
								cntr = 0;
								np_prev = np;
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
				if(thread==0) pb.checkInterrupt(); 
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
	
	
	// define gc_scoretab
	ivec gc_scoretab = ivec2_iter(gscoretab);	

	
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
					std::shuffle(itr + gc_scoretab[i], itr + gc_scoretab[i+1], lrng);
				}		
			}
		}
	}
	
	return Rcpp::List::create(Named("theta") = theta, Named("prior_log") = prior_log);
}
