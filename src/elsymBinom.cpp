#include <RcppArmadillo.h>
#include "elsymBinom.h"

using namespace arma;
using Rcpp::Named;

/*******************************************************************************
* Copies of elsym for use with lbinomial
*
* range of lbinom has to be at least three more than the possible max score
*
********************************************************************************/


void b_elsym_i_dich(const arma::mat& lbinom,const arma::ivec& a, const arma::vec& b, arma::mat& g, arma::vec& g_all, 
					const int from_item=0, const bool compute_full_elsym=true)
{
	const int nit = b.n_elem;
	int mi,sa;

	double bi;

	ivec max_score(nit, fill::zeros);
	g.zeros();
	g.at(0,0)=1;
	
	vec gw(g.n_rows);
	
	int m=0;
	for(int i=0; i<nit-1;i++)
	{
		max_score[i] = m;
		mi = m+a[i];
		
		for(int s=0; s<=m; s++)
		{
			g.at(s,i+1) += g.at(s,i) * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
			g.at(s+a[i],i+1) += g.at(s,i) * b[i] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a[i]));		
		}
		
		m = mi;
	}	
	
	// full elsym
	if(compute_full_elsym)
	{
		g_all.zeros();
		const int i = nit-1;
		mi = m+a[i];
		for(int s=0; s<=m; s++)
		{
			g_all[s] += g.at(s,i) * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
			g_all[s+a[i]] += g.at(s,i) * b[i] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a[i]));		
		}
	}
	
	int m4;
	for(int i=nit-1;i>from_item;i--)
	{
		bi=b[i];
		for(int j=0;j<i;j++) 
		{
			gw = g.col(j);
			g.col(j).zeros();
			
			m = max_score[j];
			mi = m + a[i];
			m4 = m - (m%4);
			sa = a[i];

			for(int s=0; s<=m4; s+=4, sa+=4)
			{
				g.at(s,j) += gw[s] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
				g.at(s+1,j) += gw[s+1] * std::exp(lbinom.at(m,s+1) - lbinom.at(mi,s+1));
				g.at(s+2,j) += gw[s+2] * std::exp(lbinom.at(m,s+2) - lbinom.at(mi,s+2));
				g.at(s+3,j) += gw[s+3] * std::exp(lbinom.at(m,s+3) - lbinom.at(mi,s+3));
				
				g.at(sa,j) += gw[s] * bi * std::exp(lbinom.at(m,s) - lbinom.at(mi,sa));	
				g.at(sa+1,j) += gw[s+1] * bi * std::exp(lbinom.at(m,s+1) - lbinom.at(mi,sa+1));	
				g.at(sa+2,j) += gw[s+2] * bi * std::exp(lbinom.at(m,s+2) - lbinom.at(mi,sa+2));	
				g.at(sa+3,j) += gw[s+3] * bi * std::exp(lbinom.at(m,s+3) - lbinom.at(mi,sa+3));	
	
			}
		}
		max_score += a[i];
	}
}



void b_elsym_i_poly(const arma::mat& lbinom,const arma::ivec& a,const arma::vec& b,  
					const arma::ivec& first, const arma::ivec& last,
					arma::mat& g, arma::vec& g_all,
					const int from_item=0, const bool compute_full_elsym=true)
{
	const int nit = first.n_elem;
	int mi,a1;
	double b1;
	vec gw(g.n_rows);
	
	ivec max_score(nit, fill::zeros);
	
	g.zeros();
	g.row(0).ones();	
	
	// setup	
	// partial
	int m=0;
	for(int i=0; i<nit-1;i++)
	{
		max_score[i] = m;
		mi = m + a[last[i]];

		a1 = a[first[i]];
		b1 = b[first[i]];
		
		for(int s=m;s>=0;s--)
		{				
			g.at(s,i+1)  = g.at(s,i)  * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
			g.at(s+a1,i+1) += g.at(s,i) * b1 * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a1));			
		}
		
		for(int j=first[i]+1;j<=last[i]; j++)
			for(int s=0; s<=m;s++)
				g.at(s+a[j],i+1) += g.at(s,i) * b[j] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a[j]));	
		
		m+=a[last[i]];
	}	
	
	// bonus full elsym
	if(compute_full_elsym)
	{
		const int i = nit-1;
		mi = m + a[last[i]];

		a1 = a[first[i]];
		b1 = b[first[i]];
		for(int s=m;s>=0;s--)
		{				
			g_all[s]  = g.at(s,i)  * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
			g_all[s+a1] += g.at(s,i) * b1 * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a1));			
		}
		
		for(int j=first[i]+1;j<=last[i]; j++)
			for(int s=0; s<=m;s++)
				g_all[s+a[j]] += g.at(s,i) * b[j] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a[j]));			
	}
	
	// elsym minus i
	int m4;
	for(int i=nit-1;i>from_item;i--)
	{
		
		b1 = b[first[i]];

		// for all columns
		for(int j=0;j<i;j++) 
		{
			m = max_score[j];
			m4 = m-m%4;
			mi = m+a[last[i]];
			gw = g.col(j);	
			g.col(j).zeros();
			
			
			for(int s=0, sa=a[first[i]]; s<=m4; s+=4, sa+=4)
			{
				g.at(s,j) += gw[s] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
				g.at(s+1,j) += gw[s+1] * std::exp(lbinom.at(m,s+1) - lbinom.at(mi,s+1));
				g.at(s+2,j) += gw[s+2] * std::exp(lbinom.at(m,s+2) - lbinom.at(mi,s+2));
				g.at(s+3,j) += gw[s+3] * std::exp(lbinom.at(m,s+3) - lbinom.at(mi,s+3));
				
				g.at(sa,j) += gw[s] * b1 * std::exp(lbinom.at(m,s) - lbinom.at(mi,sa));	
				g.at(sa+1,j) += gw[s+1] * b1 * std::exp(lbinom.at(m,s+1) - lbinom.at(mi,sa+1));	
				g.at(sa+2,j) += gw[s+2] * b1 * std::exp(lbinom.at(m,s+2) - lbinom.at(mi,sa+2));	
				g.at(sa+3,j) += gw[s+3] * b1 * std::exp(lbinom.at(m,s+3) - lbinom.at(mi,sa+3));	
			}
			
			for(int ii=first[i]+1; ii<=last[i]; ii++)
			{
				for(int s=0,sa=a[ii]; s<=m4; s+=4, sa+=4)
				{
					g.at(sa,j) += gw[s] * b[ii] * std::exp(lbinom.at(m,s) - lbinom.at(mi,sa));	
					g.at(sa+1,j) += gw[s+1] * b[ii] * std::exp(lbinom.at(m,s+1) - lbinom.at(mi,sa+1));	
					g.at(sa+2,j) += gw[s+2] * b[ii] * std::exp(lbinom.at(m,s+2) - lbinom.at(mi,sa+2));	
					g.at(sa+3,j) += gw[s+3] * b[ii] * std::exp(lbinom.at(m,s+3) - lbinom.at(mi,sa+3));	
				}
			}
		}
		max_score += a[last[i]];
	}
}


// bookkeeps and redirects

void elsym_i_binom(const arma::mat& lbinom, const arma::vec& b_in, const arma::ivec& a_in, 
			int* ptr_first, int* ptr_last, 
			const int nit_all, const int exclude_item,
			arma::mat& g, arma::vec& g_all, // output: g
			const int from_item, const bool compute_full_elsym)		
{

	// bookkeeping
	//const int nit_all = first_in.n_elem;
	const int nit = nit_all - (int)(exclude_item>=0);
	
	int npars=0;
	for(int i=0; i<nit_all; i++) if(i!=exclude_item)
	{
		npars +=  *(ptr_last+i) - *(ptr_first+i) + 1;
	}
	
	ivec first(nit), last(nit), a(npars);
	vec b(npars);
	
	const bool dichotomous = (npars == nit);
	
	// all contiguous
	int new_pos=0, new_i=0;
	for(int i=0;i<nit_all;i++) if(i!=exclude_item)
	{
		first[new_i] = new_pos;
		for(int j=*(ptr_first+i); j<=*(ptr_last+i); j++)
		{
			b[new_pos] = b_in[j];
			a[new_pos++] = a_in[j];
		}
		last[new_i++] = new_pos-1;
	}
	
	if(dichotomous)
		b_elsym_i_dich(lbinom, a, b, g, g_all, from_item, compute_full_elsym);
	else
		b_elsym_i_poly(lbinom, a, b, first, last, g, g_all, from_item, compute_full_elsym);	
}





void elsym_binom(const arma::mat& lbinom, const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, arma::vec& g, arma::vec& gw, const int omit_item)
{
	g.zeros();
	g[0] = 1;
	
	int m=0,mi,a1;
	double b1;
	for(int i=0; i<nit;i++) if(i != omit_item)
	{
		mi = m + a[*(last+i)];
		gw.subvec(0,m+1) = g.subvec(0,m+1);
		g.subvec(0,m+1).zeros();
		a1 = a[*(first+i)];
		b1 = b[*(first+i)];
		for(int s=0;s<=m;s++)
		{				
			g[s] += gw[s] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s));
			g[s+a1] += gw[s] * b1 * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a1));			
		}
		
		for(int j=1+*(first+i); j<=*(last+i); j++)
			for(int s=0;s<=m;s++)
				g[s+a[j]] += gw[s] * b[j] * std::exp(lbinom.at(m,s) - lbinom.at(mi,s+a[j]));			

		m = mi;
	}
}



//[[Rcpp::export]]
arma::vec elsym_binomC(const arma::mat& lbinom, const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last, const int omit_item=-1)
{
	const int nit = first.n_elem;
	
	int n_scores=1;
	for(int i=0;i<nit;i++)
		n_scores += a[last[i]];
	
	vec g(n_scores), gw(n_scores);

	elsym_binom(lbinom, b, a, first.memptr(), last.memptr(), nit, g, gw,omit_item);
	
	return g;
}



//for testing
//[[Rcpp::export]]
Rcpp::List list_elsymi_binomC(const arma::mat& lbinom, const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last)
{
	const int nit = first.n_elem;
	int ms=0;
	for(int i=0;i<nit;i++)
		ms+=a[last[i]];
	
	mat gi(ms+4,nit);
	vec g(ms+4);
	
	elsym_i_binom(lbinom,b,a,first.memptr(),last.memptr(),nit,-1,gi,g);
	
	return Rcpp::List::create(Named("gi")=gi, Named("gfull")=g);

}




