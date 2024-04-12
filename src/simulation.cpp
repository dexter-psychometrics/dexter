#include <RcppArmadillo.h>
#include "ability.h"

using namespace arma;


// [[Rcpp::export]]
arma::ivec sampleNRM_testC(	const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
							const arma::ivec& first, const arma::ivec& last)
{
	const int nit = first.n_elem;
	const int m = theta.n_elem;  
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	double u;
	int k;
	
	vec p(maxA+3);  
	vec lookup(maxA+1);
  
	ivec score(m, fill::zeros);  
  
	lookup[0] = 1;
	p[0]=1;

	for (int pers=0;pers<m;pers++)
	{
		for(int i=1;i<=maxA;i++){ lookup[i] = std::exp(i*theta[pers]);}
		for (int i=0;i<nit;i++)
		{
			k=1;
			for (int j=first[i];j<=last[i];j++) 
			{
				p[k]=p[k-1]+b[j]*lookup[a[j]];
				k++;
			}
			u=p[k-1]*R::runif(0,1);
			k=0;
			while (u>p[k]) {k++;}
			if (k>0) score[pers] += a[first[i]+k-1];
		}
	}

	return score;
}




// [[Rcpp::export]]
arma::imat sampleNRM_itemC(const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
						  const arma::ivec& first, const arma::ivec& last)
{

	const int nit = first.n_elem;
	const int m = theta.n_elem;  
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	int k;
	double u;
	
	vec p(maxA+3);
	vec lookup(maxA+1);
  
	imat score(m, nit);  
  
	lookup[0] = 1;
	p[0] = 1;

	for (int pers=0;pers<m;pers++)
	{
		for(int i=1; i<=maxA; i++)
			lookup[i] = exp(i*theta[pers]);
    
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
			if(k>0) score.at(pers, i) = a[first[i]+k-1];
		}
	}

  return score;
}


// TO DO: not tested IM yet

// [[Rcpp::export]]
arma::imat sampleIMC(const arma::vec& bIM, const arma::vec& cIM, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last,
					const arma::ivec& scoretab)
{
	// scoretab must already be based on the sample and must include 0 scores even for impossible
	const int nit = first.n_elem;
	const int max_score = scoretab.n_elem - 1;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int nP = accu(scoretab);
	const int maxsec = 200;
	const double acc = 1e-8;
	
	double u, xl, fl, f, dx, theta = -2;
	int score, k, pi;
	
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
		
		// derive theta
		xl = theta + 0.5;
		fl = Escore_single(xl, b, a, first, last, nit, maxA);
		f = Escore_single(theta, b, a, first, last, nit, maxA);
		// using secant
		for(int iter=0; iter<maxsec; iter++)
		{
			dx = (xl-theta) * (f-s)/(f-fl);
			xl = theta;
			fl = f;
			theta += dx;
			f = Escore_single(theta, b, a, first, last, nit, maxA);
			if(std::abs(dx) < acc)
				break;
		}
		
		for(int i=1; i<=maxA; i++)
			lookup[i] = exp(i*theta);
				
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
	}

	return out;
}
