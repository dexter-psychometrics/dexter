#ifndef DEXTER_PRIORS_
#define DEXTER_PRIORS_

#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <dqrng_distribution.h>

double rgamma(dqrng::xoshiro256plus& lrng, const double alpha, const double ibeta);

class hnorm_prior
{
private:
	dqrng::uniform_distribution runif;
	dqrng::normal_distribution rnorm;
	
	int npop;
	double mu, J, tau2;

public:
	arma::vec theta;
	double sigma;
	
	hnorm_prior(const arma::vec& theta_start, const double sigma_start);
	
	arma::vec as_vec();

	void update(dqrng::xoshiro256plus& lrng, const arma::vec& pv, 
					 const arma::ivec& scoretab_pop, const arma::ivec& scoretab_np, const arma::ivec& scoretab_cnp);

};


class mixture_prior
{
private:
	dqrng::uniform_distribution runif;
	dqrng::normal_distribution rnorm;
	
	double p;
public:
	arma::vec mu, sigma;
	arma::ivec pop;

	mixture_prior(dqrng::xoshiro256plus& lrng, const double p_start, const arma::vec& mu_start, const arma::vec& sigma_start, const double n);
	
	arma::vec as_vec();

	void upd_normal(dqrng::xoshiro256plus& lrng, const arma::vec& pv);

	void upd_jeffreys(dqrng::xoshiro256plus& lrng, const arma::vec& pv);

};

#endif