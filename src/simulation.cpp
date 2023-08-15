#include <RcppArmadillo.h>
#include "ability.h"

using namespace arma;

#define SEED std::round(R::runif(0,1) * 2147483647)

// cat 0 must be present
// [[Rcpp::export]]
arma::ivec sampleNRM2_test(	const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
							const arma::ivec& first, const arma::ivec& last)
{

  double u;
  const int nI = first.n_elem;
  const int m = theta.n_elem;  
  const int maxA = max(a(conv_to<uvec>::from(last))); 
  int k=0;
  vec p(maxA+3);
  vec lookup(maxA+1);
  
  ivec score(m, fill::zeros);
  
  
  lookup[0] = 1.0;

  for (int pers=0;pers<m;pers++)
  {
	for(int i=1;i<=maxA;i++){ lookup[i] = exp(i*theta[pers]);}
    for (int i=0;i<nI;i++)
    {
      p[0]=b[first[i]]; //*exp(a[first[i]]*theta[pers]); 
      k=1;
      for (int j=first[i]+1;j<=last[i];j++) // note the +1
      {
        p[k]=p[k-1]+b[j]*lookup[a[j]];
        k++;
      }
      u=p[k-1]*R::runif(0,1);
      k=0;
      while (u>p[k]) {k++;}
      if (k>0) {score[pers]+=a[first[i]+k];}
    }
  }

  return score;
}


// [[Rcpp::export]]
arma::imat sampleNRM2_item(	const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
							const arma::ivec& first, const arma::ivec& last)
{

  double u;
  const int nI = first.n_elem;
  const int m = theta.n_elem;  
  const int maxA = max(a(conv_to<uvec>::from(last))); 
  int k=0;
  vec p(maxA+3);
  vec lookup(maxA+1);
  
  imat score(m, nI);  
  
  lookup[0] = 1.0;

  for (int pers=0;pers<m;pers++)
  {
	for(int i=1; i<=maxA; i++)
		lookup[i] = exp(i*theta[pers]);
    
	for (int i=0;i<nI;i++)
    {
		p[0] = b[first[i]]; 
		k = 1;
		for (int j=first[i]+1; j<=last[i]; j++) // note the +1
		{
			p[k] = p[k-1]+b[j]*lookup[a[j]];
			k++;
		}
		u = p[k-1]*R::runif(0,1);
		k = 0;
		while(u>p[k]) k++;
		score.at(pers, i) = a[first[i]+k];
    }
  }

  return score;
}


// [[Rcpp::export]]
arma::imat sampleIM(const arma::vec& bIM, const arma::vec& cIM, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last,
					const arma::ivec& scoretab)
{
	// scoretab must already be based on the sample and must include 0 scores even for impossible
	const int nI = first.n_elem;
	const int max_score = scoretab.n_elem - 1;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int nP = accu(scoretab);
	const int maxsec = 200;
	const double acc = 1e-8;
	double theta = -2;
	
	vec logb = log(bIM), logc = log(cIM);
	
	vec b(&bIM[0], bIM.n_elem);

	ivec cs_scoretab(max_score+1);
	vec lookup(maxA+1);
	vec p(maxA+3);
	
	imat out(nP, nI);	
	
	lookup[0] = 1;
	
	cs_scoretab[0] = 0;
	for(int i=0; i<max_score; i++)
		cs_scoretab[i+1] = scoretab[i] + cs_scoretab[i];

	
	for(int s=1; s<max_score; s++)
	{
		for(int i=0; i<nI; i++)
			for (int j=first[i]+1; j<=last[i]; j++) 
				b[j] = exp(logb[j] + s * a[j] * logc[i]);
		
		// derive theta
		double xl = theta + 0.5;
		double fl = Escore_single(xl, b, a, first, last, nI, maxA),
				f = Escore_single(theta, b, a, first, last, nI, maxA);
		// using secant
		for(int iter=0; iter<maxsec; iter++)
		{
			double dx = (xl-theta) * (f-s)/(f-fl);
			xl = theta;
			fl = f;
			theta += dx;
			f = Escore_single(theta, b, a, first, last, nI, maxA);
			if(std::abs(dx) < acc)
				break;
		}
		
		for(int i=1; i<=maxA; i++)
			lookup[i] = exp(i*theta);
				
		int pi = cs_scoretab[s];
		// sample
		while(pi < cs_scoretab[s+1])
		{
			int score = 0;
			for (int i=0;i<nI;i++)
			{				
				p[0] = b[first[i]]; 
				int k = 1;				
				for (int j=first[i]+1; j<=last[i]; j++) 
				{
					p[k] = p[k-1]+b[j]*lookup[a[j]];
					k++;
				}
				double u = p[k-1]*R::runif(0,1);
				k = 0;
				while(u>p[k]) k++;
				score += a[first[i]+k];
				out.at(pi, i) = a[first[i]+k];
				if(score > s )
					break;
			}
			if(score == s)
				pi++;	
		}
	}

	return out;
}
