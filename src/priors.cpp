#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>
#include "priors.h"
#include "shared.h"

#define SEED std::round(R::runif(0,1) * 2147483647)

using namespace arma;

/* *********************************************************************************************************
* Random deviate for gamma not available in dq
* 
* To do: check if available in boost
*
************************************************************************************************************/

double rgamma(dqrng::xoshiro256plus& lrng, const double alpha, const double ibeta)
{
	const double beta = 1/ibeta;
	const double alph = alpha + (alpha < 1);
	const double a1 = alph - 1.0/3;
	const double a2 = 1/std::sqrt(9*a1);
	
	dqrng::uniform_distribution runif;
	dqrng::normal_distribution rnorm;
		
	// marsaglia & Tsang, via numerical recipes
	double u,v,x;
	do
	{
		do
		{
			x = rnorm(lrng);
			v = 1 + a2*x;
			
		} while (v<=0);			
			
		v = CUB(v);
		u = runif(lrng);
			
	} while(u > 1-0.331 * QRT(x) && std::log(u) > 0.5 * SQR(x) + a1*(1-v+std::log(v)));
			
	if(alph > alpha)
	{
		u = runif(lrng);
		while(u==0)
			u = runif(lrng);
		return std::pow(u, 1/alpha)*a1*v/beta;
	}
	else
	{
		return a1*v/beta;
	}	
}
	

// transformations of rgamma
	
double rchisq(dqrng::xoshiro256plus& lrng, const double df)
{
	return 2 * rgamma(lrng, df/2, 1.0);
}
	
double rinvchisq(dqrng::xoshiro256plus& lrng, const double df, const double scale)
{
	double z = rchisq(lrng, df);
	if(z==0) z= 1e-100;
	return (df*scale) / z;
}

double rbeta(dqrng::xoshiro256plus& lrng, const double alpha, const double beta)
{
	double x = rgamma(lrng, alpha, 1.0);
	double y = rgamma(lrng, beta, 1.0);
	
	return x/(x+y);
}

// for testing

// [[Rcpp::export]]
arma::vec test_rgamma(const int n, const double alpha, const double beta)
{
	vec out(n);
	dqrng::xoshiro256plus rng(SEED);
	for(int i=0; i<n; i++)
	{
		out[i] = rgamma(rng, alpha, beta);	
	}
	return out;
}


// [[Rcpp::export]]
arma::vec test_rinvchisq(const int n, const double df, const double scale)
{
	vec out(n);
	dqrng::xoshiro256plus rng(SEED);
	for(int i=0; i<n; i++)
	{
		out[i] = rinvchisq(rng, df, scale);	
	}
	return out;
}

// [[Rcpp::export]]
arma::vec test_rbeta(const int n, const double alpha, const double beta)
{
	vec out(n);
	dqrng::xoshiro256plus rng(SEED);
	for(int i=0; i<n; i++)
	{
		out[i] = rbeta(rng, alpha, beta);	
	}
	return out;
}



/* *********************************************************************************************************
* Prior updating functions for plausible values that can be used in parallel
* 
* 
*
************************************************************************************************************/

// constructor
hnorm_prior::hnorm_prior(const vec& theta_start, const double sigma_start)
{
	npop = theta_start.n_elem;
	J = (double)npop;
	theta = theta_start;
	mu = mean(theta);
	sigma = sigma_start;
	runif = dqrng::uniform_distribution(0, 1);
	rnorm = dqrng::normal_distribution(0, 1);
}	
	

// The variable names in Gelman, Rubin, et. al are used here
// therefore theta are the group means and sigma is the common group standard deviation
void hnorm_prior::update(dqrng::xoshiro256plus& lrng, const vec& pv, const ivec& scoretab_pop, const ivec& scoretab_np, const ivec& scoretab_cnp)
{
	const double N = (double)pv.n_elem;
		
	if(npop==1)
	{
		// Single normal prior in case of one group
		const double m = mean(pv);
		const double v = var(pv-m);
		sigma = std::sqrt(1/rgamma(lrng, (N-1)/2, 1/(((N-1)/2)*v)));
		theta[0] = rnorm(lrng) * (sigma/std::sqrt(N)) + m;	
	}
	else
	{
		// data summary, walk over pv's per scoretab to produce y_bar and sigma_hat2
		
		const int ntab = scoretab_pop.n_elem;
		
		int pop;			
		double theta_hat, V_theta, mu_hat, tau_hat2, sigma2, sigma_hat2 = 0;
		
		ivec n(npop, fill::zeros);
		vec ybar(npop,fill::zeros);

		for(int tab=0;tab<ntab;tab++)
		{
			pop = scoretab_pop[tab];
			n[pop] += scoretab_np[tab];
			for(int p=scoretab_cnp[tab]; p<scoretab_cnp[tab+1]; p++)
			{
				ybar[pop] += pv[p];
				sigma_hat2 += SQR(pv[p] - theta[pop]);
			}
		}
		sigma_hat2 /= N;
		ybar /= conv_to<vec>::from(n);
		
		// The sampling		
		if(npop<5)
		{
			tau2 = 1/rgamma(lrng,1.0,1.0);
		}
		else
		{		
			tau_hat2 = accu(square(theta-mu))/(J-1);
			tau2 = rinvchisq(lrng, J-1, tau_hat2);
		} 


		sigma2 = rinvchisq(lrng, N, sigma_hat2);
		
		mu_hat = mean(theta);
		mu = rnorm(lrng) * std::sqrt(tau2/J) + mu_hat;
		
		for(int j=0; j<npop; j++)
		{
			V_theta = 1/(1/tau2 + n[j]/sigma2); 
			
			theta_hat = (mu/tau2 + ybar[j]*n[j]/sigma2)/(1/tau2+n[j]/sigma2);
			theta[j] = rnorm(lrng) * std::sqrt(V_theta) + theta_hat;
		}		
		sigma = std::sqrt(sigma2);
	}
}

// all priors as one vector for inclusion in log
vec hnorm_prior::as_vec()
{
	vec v(theta.n_elem + 3);
	v.tail(theta.n_elem) = theta;
	v[0] = mu;
	v[1] = sigma;
	v[2] = std::sqrt(tau2);
	return v;
}

// Mixture prior


mixture_prior::mixture_prior(dqrng::xoshiro256plus& lrng, const double p_start, const vec& mu_start, const vec& sigma_start, const double n)
{
	runif = dqrng::uniform_distribution(0, 1);
	rnorm = dqrng::normal_distribution(0, 1);
	p = p_start;
	mu = vec(2);
	sigma = vec(2);
	mu = mu_start;
	sigma = sigma_start;
	
	pop = ivec(n);
	for(int i=0;i<n;i++)
		pop[i] = runif(lrng) > 0.5;
}

vec mixture_prior::as_vec()
{
	vec out = {p, mu[0], mu[1], sigma[0], sigma[1]};
	return out;	
}
  
  
// to do: lit nakijken op typfout
void mixture_prior::upd_normal(dqrng::xoshiro256plus& lrng, const vec& pv)
{
	const int n = pv.n_elem;
	double p1, p2;
	const double mean_pv = mean(pv), var_pv = var(pv);
	
	ivec popn(2,fill::zeros);
	vec popmean(2,fill::zeros), popvar(2,fill::zeros);
	
	// hyper prior
	const double L=1;
	const double V=3;
	
	// latent group membership
	for(int i=0; i<n; i++)
	{
		p1 = p * R::dnorm(pv[i], mu[0], sigma[0], false);
		p2 = (1-p) * R::dnorm(pv[i], mu[1], sigma[1], false);
		
		if(p1>0) pop[i] = runif(lrng) > p1/(p1+p2);  // > larger than, otherwise p keeps switching
		else pop[i] = runif(lrng) > 0.5;	
		
		popn[pop[i]]++;
		popmean[pop[i]] += pv[i];
	}
	popmean /= conv_to<vec>::from(popn);
	
	// means
	for(int j=0; j<2; j++)
	{
		double m = (L * mean_pv + popn[j] * popmean[j]) / (L + popn[j]);
		double s = std::sqrt(SQR(sigma[j]) / (L + popn[j]));
		mu[j] = rnorm(lrng) * s + m;
	}
	
	// variances
	for(int i=0; i<n; i++) popvar[pop[i]] += SQR(pv[i] - popmean[pop[i]]);
	popvar /= conv_to<vec>::from(popn-1);
	
	for(int j=0; j<2; j++) if(popn[j] > 1)
	{
		double shape = 0.5 * (V + popn[j]); 
		double rate = 0.5 * (var_pv + popn[j] * popvar[j] +  (L * popn[j]/(L + popn[j])) * SQR(mean_pv - popmean[j]));
		
		sigma[j] = std::sqrt(1/rgamma(lrng, shape, 1/rate));	
	}
	
    p = rbeta(lrng, popn[0] + 2.0, popn[1] + 2.0);
}


void mixture_prior::upd_jeffreys(dqrng::xoshiro256plus& lrng, const vec& pv)
{
	const int n = pv.n_elem;
	double p1, p2;
	
	ivec popn(2,fill::zeros);
	vec popmean(2,fill::zeros), popsumsq(2,fill::zeros);
	
	// latent group membership
	for(int i=0; i<n; i++)
	{
		p1 = p * R::dnorm(pv[i], mu[0], sigma[0], false);		
		p2 = (1-p) * R::dnorm(pv[i], mu[1], sigma[1], false);
		
		if(p1>0) pop[i] = runif(lrng) <= p1/(p1+p2);
		else pop[i] = runif(lrng) > 0.5;	
		
		popn[pop[i]]++;
		popmean[pop[i]] += pv[i];
	}
	popmean /= conv_to<vec>::from(popn-1);
	
	// means
	for(int j=0; j<2; j++)
		mu[j] = rnorm(lrng) * sigma[j]/std::sqrt(popn[j]) + popmean[j];
	
	
	// variances
	for(int i=0; i<n; i++) popsumsq[pop[i]] += SQR(pv[i] - popmean[pop[i]]);

	for(int j=0; j<2; j++) if(popn[j] > 1)
		sigma[j] = 1/std::sqrt(rgamma(lrng, popn[j]/2, 1/(popsumsq[j])/2));
		
	
	// membership probability
	p = rbeta(lrng, popn[0] + 0.5, popn[1] + 0.5);
}


