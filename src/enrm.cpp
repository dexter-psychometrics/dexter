
#include <RcppArmadillo.h>
#include "elsym.h"
#include "elsymBinom.h"
#include "shared.h"
#include "myomp.h"

#define ELSYM_USE_LD false

#if ELSYM_USE_LD
#define ELS_VEC ldvec
#define ELS_MAT ldmat
#else
#define ELS_VEC vec
#define ELS_MAT mat
#endif

using namespace arma;

/*********************************************************************
* CALIBRATION
*
* using elsym_i
* everything is without the 0 category
*
**********************************************************************/

//[[Rcpp::export]]
arma::vec Expect(const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit)
{
	const int nbooklets = nit.n_elem;
	const int npar = b.n_elem;
	const int max_nit = max(nit);
	const int max_nscore = max(n_score);	
	
	ELS_MAT gi(max_nscore+3,max_nit);
	ELS_VEC g_all(max_nscore);
	
	vec E(npar, fill::zeros);
	
	// cumulatives for bookkeeping
	ivec cnit(nbooklets), cn_score(nbooklets); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nbooklets-1) = cumsum(nit.head(nbooklets-1));
	cn_score.tail(nbooklets-1) = cumsum(n_score.head(nbooklets-1));	
	
	
	for(int bk=0; bk < nbooklets; bk++)
	{
		int ipos = cnit[bk];
		int spos = cn_score[bk];
		
		elsym_i(b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], -1, gi, g_all);		
		
		
		for(int i=0; i<nit[bk]; i++, ipos++)
		{
			for(int j=first[ipos]; j<=last[ipos]; j++)
			{
				for(int s = a[j]; s < n_score[bk]; s++) if(g_all[s]>0)
				{
					E[j] += scoretab[spos + s] * b[j] * (gi.at(s-a[j],i)/ g_all[s]);
				}
			}
		}	
	}
	return E;
}




//[[Rcpp::export]]
arma::vec Expect_binom(const arma::mat& lbinom, const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit)
{
	const int nbooklets = nit.n_elem;
	const int npar = b.n_elem;
	const int max_nit = max(nit);
	const int max_nscore = max(n_score);
	int mx1,mx2;
	
	mat gi(max_nscore+3,max_nit);
	vec g_all(max_nscore);
	
	vec E(npar, fill::zeros);
	
	// cumulatives for bookkeeping
	ivec cnit(nbooklets), cn_score(nbooklets); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nbooklets-1) = cumsum(nit.head(nbooklets-1));
	cn_score.tail(nbooklets-1) = cumsum(n_score.head(nbooklets-1));	
	
	
	for(int bk=0; bk < nbooklets; bk++)
	{
		int ipos = cnit[bk];
		int spos = cn_score[bk];
		
		elsym_i_binom(lbinom,b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], -1, gi, g_all);		
		
		mx1 = n_score[bk]-1;
		
		for(int i=0; i<nit[bk]; i++, ipos++)
		{
			mx2 = mx1 - a[last[ipos]];
			for(int j=first[ipos]; j<=last[ipos]; j++)
			{
				for(int s = a[j]; s < n_score[bk]; s++) if(g_all[s]>0)
					E[j] += scoretab[spos + s] * b[j] * (gi.at(s-a[j],i)/ g_all[s]) * std::exp(lbinom.at(mx2,s-a[j]) - lbinom.at(mx1,s)); 
			}
		}	
	}
	return E;
}



//[[Rcpp::export]]
void Hess(const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit, const int max_cores,
				 arma::vec& E, arma::mat& H)
{
	const int nbooklets = nit.n_elem;
	const int npar = b.n_elem;
	const int max_nit = max(nit);
	const int max_nscore = max(n_score);
	
	E.zeros();
	H.zeros();
	
	// cumulatives for bookkeeping
	ivec cnit(nbooklets), cn_score(nbooklets); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nbooklets-1) = cumsum(nit.head(nbooklets-1));
	cn_score.tail(nbooklets-1) = cumsum(n_score.head(nbooklets-1));	
	
	
	

#pragma omp parallel num_threads(max_cores)
{
	ELS_MAT gi(max_nscore+3,max_nit);
	ELS_MAT gik(max_nscore+3,max_nit);
	
	ELS_VEC g_all(max_nscore);
	ELS_VEC dummy_g(1);
	
	vec d1_part(max_nscore);
	
	double tmp;
	int ipos1;

#pragma omp for reduction(+: E,H)
	for(int bk=0; bk < nbooklets; bk++)
	{
		ipos1 = cnit[bk];
		const int spos = cn_score[bk];
		const int bmax = n_score[bk]-1;
		
		elsym_i(b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], -1, gi, g_all, 0, true);
		
		// item1 is positie item i binnen boekje en binnen gi		
		// ipos1 en ipos2 worden gebruikt om first en last te indexeren
		
		for(int item1 = 0; item1<nit[bk]; item1++,ipos1++)
		{
			elsym_i(b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], item1, gik, dummy_g, item1, false);
			
			for(int j=first[ipos1]; j<=last[ipos1]; j++)
			{
				// diagonal
				for (int s=a[j]; s<=bmax; s++) if(g_all[s]>0)
				{
					tmp = b[j] * ( gi.at(s-a[j],item1)/ g_all[s] );
					d1_part[s] = scoretab[s+spos ] * tmp;		
					E[j] += d1_part[s]; 
					H.at(j,j) += d1_part[s] * (1-tmp);				
				}
				// block diagonal, poytomous
				for (int k=(j+1);k<=last[ipos1];k++)
				{
					for (int s=a[k]; s<=bmax; s++) if(g_all[s]>0)
					{
						H.at(j,k) -= d1_part[s] * b[k] * (gi.at(s-a[k],item1)/ g_all[s]); 
					}				
				}
				// item2 is positie item k binnen gik
				for(int item2 = item1+1, ipos2=ipos1+1; item2<nit(bk); item2++, ipos2++)
				{
					for(int k = first[ipos2]; k <= last[ipos2]; k++)
					{
						for (int s=a[j];s<=bmax;s++) if(g_all[s]>0)
						{							
							if (s >= (a[j]+a[k]))
							{
								H.at(j,k) += scoretab[s+spos] * b[j] * b[k] * (gik.at(s-a[j]-a[k],item2-1)/ g_all[s]);
							}
							if(s>=a[j] && s>=a[k])
							{
								H.at(j,k) -= d1_part[s] * b[k] * (gi.at(s-a[k],item2)/ g_all[s]); 
							}
						}
					}
				}			
			}		
		}	
	}
}	
	for(int i=0;i<npar;i++)
		for(int j=i+1; j<npar;j++)
			H.at(j,i) = H.at(i,j);
	
}



//[[Rcpp::export]]
void Hess_binom(const arma::mat& lbinom,const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit, const int max_cores,
				 arma::vec& E, arma::mat& H)
{
	const int nbooklets = nit.n_elem;
	const int npar = b.n_elem;
	const int max_nit = max(nit);
	const int max_nscore = max(n_score);
	
	E.zeros();
	H.zeros();
	
	// cumulatives for bookkeeping
	ivec cnit(nbooklets), cn_score(nbooklets); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nbooklets-1) = cumsum(nit.head(nbooklets-1));
	cn_score.tail(nbooklets-1) = cumsum(n_score.head(nbooklets-1));	
		

#pragma omp parallel num_threads(max_cores)
{
	mat gi(max_nscore+3,max_nit);
	mat gik(max_nscore+3,max_nit);
	
	vec g_all(max_nscore);
	vec dummy_g(1);
	
	vec d1_part(max_nscore);
	
	double tmp;	
	int mx1, mx2, mx3, ipos1;

#pragma omp for reduction(+: E,H)
	for(int bk=0; bk < nbooklets; bk++)
	{
		ipos1 = cnit[bk];
		const int spos = cn_score[bk];
		const int bmax = n_score[bk]-1;
		
		elsym_i_binom(lbinom,b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], -1, gi, g_all, 0, true);
		
		// item1 is positie item i binnen boekje en binnen gi		
		// ipos1 en ipos2 worden gebruikt om first en last te indexeren
		
		for(int item1 = 0; item1<nit[bk]; item1++,ipos1++)
		{
			elsym_i_binom(lbinom,b,a, first.memptr() + cnit[bk], last.memptr() + cnit[bk], nit[bk], item1, gik, dummy_g, item1, false);
			
			mx1 = bmax-a[last[ipos1]];
			
			for(int j=first[ipos1]; j<=last[ipos1]; j++)
			{
				// diagonal
				for (int s=a[j]; s<=bmax; s++) if(g_all[s]>0)
				{
					tmp = b[j] * ( gi.at(s-a[j],item1)/ g_all[s] ) * std::exp(lbinom.at(mx1,s-a[j]) - lbinom.at(bmax,s));
					d1_part[s] = scoretab[s+spos ] * tmp;		
					E[j] += d1_part[s]; 
					H.at(j,j) += d1_part[s] * (1-tmp);				
				}
				// block diagonal, polytomous
				for (int k=(j+1);k<=last[ipos1];k++)
				{
					for (int s=a[k]; s<=bmax; s++) if(g_all[s]>0)
					{
						H.at(j,k) -= d1_part[s] * b[k] * (gi.at(s-a[k],item1)/ g_all[s]) * std::exp(lbinom.at(mx1,s-a[k]) - lbinom.at(bmax,s));
					}				
				}
				// item2 is positie item k binnen gik
				
				for(int item2 = item1+1, ipos2=ipos1+1; item2<nit(bk); item2++, ipos2++) 
				{
					mx2 = bmax-a[last[ipos2]]-a[last[ipos1]];
					mx3 = bmax - a[last[ipos2]];
					for(int k = first[ipos2]; k <= last[ipos2]; k++)
					{
						for (int s=a[j];s<=bmax;s++) if(g_all[s]>0)
						{							
							if (s >= (a[j]+a[k]))
							{
								H.at(j,k) += scoretab[s+spos] * b[j] * b[k] * (gik.at(s-a[j]-a[k],item2-1)/ g_all[s]) * std::exp(lbinom.at(mx2,s-a[j]-a[k]) - lbinom.at(bmax,s)); 
							}
							if(s>=a[j] && s>=a[k])
							{
								H.at(j,k) -= d1_part[s] * b[k] * (gi.at(s-a[k],item2)/ g_all[s]) * std::exp(lbinom.at(mx3,s-a[k]) - lbinom.at(bmax,s)); 
							}
						}
					}
				}			
			}		
		}	
	}
}	
	for(int i=0;i<npar;i++)
		for(int j=i+1; j<npar;j++)
			H.at(j,i) = H.at(i,j);
	
}
