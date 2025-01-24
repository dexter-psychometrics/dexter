#include <RcppArmadillo.h>
#include "shared.h"
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

