
#include <RcppArmadillo.h>
#include <xoshiro.h>
#include "myomp.h"
#include "progress.h"
#include "priors.h"
#include "elsym.h"


using namespace arma;
using Rcpp::Named;

#define SEED std::round(R::runif(0,1) * 2147483647)

/*
void rdirichlet(dqrng::xoshiro256plus& lrng, const int *scoretab_ptr, double *output_ptr, const int max_score, const double prior=0.1)
{
	long double sm=0;
	
	for (int s=0; s<=max_score;s++)
	{
		*(output_ptr+s) = rgamma(lrng, prior + *(scoretab_ptr+s),1.0);
		sm += *(output_ptr+s);
	}
	for (int s=0; s<=max_score;s++)
		*(output_ptr+s) /= sm;
}
*/


// [[Rcpp::export]]
Rcpp::List calibrate_Bayes_chains(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, 
				 const arma::ivec& ib, const arma::ivec& bi, const arma::ivec& nbi, 
				 const arma::ivec& nib,
				 arma::ivec& bfirst, arma::ivec& blast,
				 const arma::ivec& bmax, const arma::ivec& m,
				 const arma::ivec& sufI, const arma::ivec& sufI_zero, const arma::ivec& bkscoretab,
				 const arma::mat& b_start,  const arma::ivec& item_fixed, 
				 const int warmup, const int step, const int ndraws,
				 const arma::ivec progress_init, const int max_cores,
				 const double prior_eta=0.5, const double prior_rho=0.5)
{
	const int nchains = b_start.n_cols;
	const int max_cat = max(last-first)+1;
	
	progress_prl pb(ndraws, progress_init);	

	dqrng::xoshiro256plus rng(SEED);
	
	// cumulatives for bookkeeping
	const int nit = nbi.n_elem;
	ivec cnbi(nit+1);
	cnbi(0) = 0;
	cnbi.tail(nit) = cumsum(nbi);
	
	const int nbk = nib.n_elem;
	ivec cnib(nbk+1);
	cnib(0) = 0;
	cnib.tail(nbk) = cumsum(nib);
	
	ivec cbmax(nbk+1);
	cbmax[0] = 0;	
	for(int k=0; k<nbk; k++)
		cbmax[k+1] = cbmax[k] + bmax[k]+1;	
	

	const int max_bscore = max(bmax);

	mat out_b(b_start.n_rows,ndraws), out_bklambda(bkscoretab.n_elem, ndraws);
	
	arma::ivec chain_start(nchains,fill::zeros);
	
#pragma omp parallel num_threads(max_cores)
	{
		const int thread = omp_get_thread_num();
		
		// working variables
		vec y(max_cat, fill::zeros), z(nbk, fill::zeros);
		vec bklambda(bkscoretab.n_elem, fill::ones);	// to do: onmogelijke scores?
		
		vec pi_k(max_bscore+1, fill::zeros), g(max_bscore+1), gw(max_bscore+1);
	
		vec b(b_start.n_rows);

		int bk,out_index;
		long double sm;
		double r_zero, b_zero;
		
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);	
		
#pragma omp for
		for(int chain=0; chain<nchains; chain++)
		{
			out_index = chain*(ndraws/nchains) + std::min(ndraws % nchains, chain);
			chain_start[chain] = out_index+1;// R index
			b = b_start.col(chain);

			const int niter = (ndraws/nchains + (chain < ndraws % nchains)-1)*step + warmup + 1;
			
			for (int iter=0; iter<niter; iter++)
			{
				if(pb.interrupted()) break;
				for (int k=0; k<nbk; k++)
				{
					elsym(b, a, bfirst.memptr() + cnib[k], blast.memptr() + cnib[k], nib[k],g,gw);

					sm = 0;
					for (int s=0; s<=bmax[k];s++)
					{
					  if (g[s]>0){
						pi_k[s] = rgamma(lrng, bkscoretab[cbmax[k]+s]+0.1, 1.0);
						sm += pi_k[s];
					  }
					}
					
					for (int s=0; s<=bmax[k];s++)
					{
					  if(g[s]>0)
					  {
						pi_k[s] = pi_k[s]/sm;
						bklambda[cbmax[k]+s] = (pi_k[s]*m[k])/g[s];
					  }
					}
					
					
					z[k] = rgamma(lrng, m[k], 1.0/m[k]);
					
				}
				// bayes_items
				// dexter lets the fixed b's vary within an iteration, weird. Also beneficial? order sensitive? should we shuffle anyways?
				for (int i=0; i<nit;i++) if(item_fixed[i] == 0) 
				{
					y.zeros();
					r_zero = 0;
					
					for (int ib_indx=cnbi[i]; ib_indx<cnbi[i+1]; ib_indx++)
					{
						bk = ib[ib_indx];
						elsym(b, a, bfirst.memptr() + cnib[bk], blast.memptr() + cnib[bk], nib[bk],g,gw, bi[ib_indx]);
							
						for (int s=0; s<=bmax[bk]-a[last[i]]; s++)
						{
							r_zero += z[bk] * g[s] * bklambda[cbmax[bk]+s];
							
							for (int j=first[i], c=0; j<=last[i]; j++,c++)
								y[c] += z[bk] * g[s] * bklambda[cbmax[bk]+s+a[j]];
						}						
					}
					b_zero = rgamma(lrng, sufI_zero[i] + prior_eta, 1/(r_zero + prior_rho));
					
					for (int j=first[i], c=0; j<=last[i]; j++,c++)
						b[j] = rgamma(lrng, sufI[j] + prior_eta, 1/(y[c] + prior_rho))/b_zero;					
					
				}
				
				b.replace(datum::nan, 1); 	

				if(iter >= warmup && (iter - warmup) % step == 0)	
				{
					pb.tick(thread == 0);
					out_b.col(out_index++) = b;
					if(thread==0) pb.checkInterrupt(); 
				}
			}
		}
	}	
	
	return Rcpp::List::create(Rcpp::Named("b") = out_b.t(), Rcpp::Named("lambda") = out_bklambda.t(), Named("chain_start") = chain_start);
}

