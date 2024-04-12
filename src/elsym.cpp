
#include <RcppArmadillo.h>
#include "shared.h"
#include "elsym.h"

using namespace arma;
using Rcpp::Named;

template <class M, class V>
void elsym_i_dich(const arma::ivec& a, const arma::vec& b, M& g, V& g_all, 
					const int from_item=0, const bool compute_full_elsym=true)
{
	const int nit = b.n_elem;
	double bi;

	ivec max_score(nit, fill::zeros);
	g.zeros();
	g.row(0).ones();	
	
	int m=0;
	for(int i=0; i<nit-1;i++)
	{
		max_score[i] = m;
		// update next column, so walk forward is allowed
		for(int s=1;s<a[i];s++) g.at(s,i+1) = g.at(s,i);
		for(int s=0; s<=m; s++)
		{
			g.at(s+a[i],i+1) = g.at(s+a[i],i) + g.at(s,i)*b[i]; // dont change to +=
		}
		
		m+=a[i];
	}	
	
	// full elsym
	if(compute_full_elsym)
	{
		g_all.zeros();
		const int i = nit-1;
		for(int s=0;s<a[i];s++) g_all[s] = g.at(s,i);
		for(int s=0;s<=m; s++)
			g_all[s+a[i]] = g.at(s+a[i],i) + g.at(s,i)*b[i];
	}
	
	int m4;
	for(int i=nit-1;i>from_item;i--)
	{
		bi = b[i];
		for(int j=0;j<i;j++) 
		{
			m = max_score[j];
			m4 = m - (m%4);

			for(int s=m4,sa=m4+a[i];s>=0;s-=4,sa-=4)
			{
				g.at(sa+3,j) += g.at(s+3,j) * bi;
				g.at(sa+2,j) += g.at(s+2,j) * bi;
				g.at(sa+1,j) += g.at(s+1,j) * bi;
				g.at(sa,j) += g.at(s,j) * bi;				
			}
		}
		max_score += a[i];
	}
}

template <class M, class V>
void elsym_i_poly(const arma::ivec& a,const arma::vec& b,  
					const arma::ivec& first, const arma::ivec& last,
					M& g, V& g_all, 
					const int from_item=0, const bool compute_full_elsym=true)
{
	const int nit = first.n_elem;
	int a1;
	double b1,bi;
	V gw(g.n_rows);
	
	ivec max_score(nit, fill::zeros);
	
	g.zeros();
	g.row(0).ones();	
	
	// setup	
	// partial
	int m=0;
	for(int i=0; i<nit-1;i++)
	{
		max_score[i] = m;

		a1 = a[first[i]];
		b1 = b[first[i]];
		for(int s=m; s>=0; s--)
		{
			g.at(s,i+1) = g.at(s,i);
			g.at(s+a1,i+1) += g.at(s,i)*b1; 
		}
		for(int j=first[i]+1;j<=last[i]; j++)
			for(int s=0; s<=m;s++)
				g.at(s+a[j],i+1) += g.at(s,i)*b[j]; 
		
		m+=a[last[i]];
	}	
	
	// bonus full elsym
	if(compute_full_elsym)
	{
		const int i = nit-1;
		g_all.zeros();
		a1 = a[first[i]];
		b1 = b[first[i]];
		for(int s=m; s>=0; s--)
		{
			g_all[s] = g.at(s,i);
			g_all[s+a1] += g.at(s,i)*b1; 
		}
		for(int j=first[i]+1;j<=last[i]; j++)
			for(int s=0; s<=m;s++)
				g_all[s+a[j]] += g.at(s,i)*b[j]; 
		
	}
	
	// elsym minus i
	for(int i=nit-1;i>from_item;i--)
	{
		// for all columns
		for(int j=0;j<i;j++) 
		{
			m = max_score[j];
			m -= m%4;
			gw = g.col(j);			
			for(int ii=first[i]; ii<=last[i]; ii++) // shortest loop
			{
				bi = b[ii];
				for(int s=0,sa=a[ii];s<=m;s+=4,sa+=4)
				{
					g.at(sa,j) += gw[s] * bi;
					g.at(sa+1,j) += gw[s+1] * bi;
					g.at(sa+2,j) += gw[s+2] * bi;
					g.at(sa+3,j) += gw[s+3] * bi;
				}
			}
		}
		max_score += a[last[i]];
	}
}




// bookkeeps and redirects
template <class M, class V>
void elsym_i(const arma::vec& b_in, const arma::ivec& a_in, 
			int* ptr_first, int* ptr_last, 
			const int nit_all, const int exclude_item,
			M& g, V& g_all, 
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
		elsym_i_dich(a, b, g, g_all, from_item, compute_full_elsym);
	else
		elsym_i_poly(a, b, first, last, g, g_all, from_item, compute_full_elsym);	
}


// vanilla simple elsym

template <class V>
void elsym(const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, V& g, V& gw, const int omit_item)
{
	g.zeros();
	g[0] = 1;
	
	int m=0;
	for(int i=0; i<nit;i++) if(i != omit_item)
	{
		//gw.subvec(0,m+1) = g.subvec(0,m+1);
		std::copy(g.memptr(),g.memptr()+m+1,gw.memptr());
		for(int j=*(first+i); j<=*(last+i); j++)
		{
			for(int s=0;s<=m;s++)
			{
				g[s+a[j]] += gw[s] * b[j];			
			}
		}
		m += a[*(last+i)];
	}
}

// explicit initiation of the long double form
template void elsym_i<ldmat,ldvec>(const arma::vec& b_in, const arma::ivec& a_in, 
			int* ptr_first, int* ptr_last, 
			const int nit_all, const int exclude_item,
			ldmat& g, ldvec& g_all, 
			const int from_item=0, const bool compute_full_elsym=true);		

// strange that this is necessary since it is clearly used in elsymC? maybe because of the default argument?
template void elsym<vec>(const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, vec& g, vec& gw, const int omit_item=-1);


//[[Rcpp::export]]
arma::vec elsymC(const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last,const int omit_item=-1)
{
	const int nit = first.n_elem;
	
	int n_scores=1;
	for(int i=0;i<nit;i++)
		n_scores += a[last[i]];
	
	vec g(n_scores,fill::zeros), gw(n_scores);

	elsym(b, a, first.memptr(), last.memptr(), nit, g, gw, omit_item);
	
	return g;
}


//for testing
//[[Rcpp::export]]
Rcpp::List list_elsymiC(const arma::vec& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last)
{
	const int nit = first.n_elem;
	int ms=0;
	for(int i=0;i<nit;i++)
		ms+=a[last[i]];
	
	mat gi(ms+4,nit);
	vec g(ms+4);
	
	elsym_i(b,a,first.memptr(),last.memptr(),nit,-1,gi,g);
	
	return Rcpp::List::create(Named("gi")=gi, Named("gfull")=g);

}





//[[Rcpp::export]]
arma::ivec possible_scores_C(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last)
{
	const int nI = last.n_elem;
	int ms = 0;
	for(int i=0; i<nI; i++) ms += a[last[i]];
	
	ivec ps(ms+1, fill::zeros);
	ps[0] = 1;
	int m = 0;
	for(int i=0; i<nI; i++)
	{
		for(int s=m; s>=0; s--)
			if(ps[s] == 1)
				for(int j=first[i]; j<=last[i]; j++)
					ps[s+a[j]] = 1;
		m += a[last[i]];
	}	
	return ps;
}


void elsym_partial(const arma::vec& b, const arma::ivec& a, int *first, int *last, const int nit, std::vector<long double>& g, std::vector<long double>& gw, const int max_score, const int omit_item=-1)
{
	std::fill(g.begin(),g.end(),0);
	g[0] = 1;
	const int mx = std::max(max_score-1,1);
	
	int m=0;
	for(int i=0; i<nit;i++) if(i != omit_item)
	{
		std::copy(g.begin(),g.begin()+m+1,gw.begin());
		for(int j=*(first+i); j<=*(last+i); j++)
		{
			for(int s=0;s<=m;s++)
			{
				g[s+a[j]] += gw[s] * b[j];			
			}
		}
		m = std::min(m+a[*(last+i)],mx);
	}
}

// Interaction model, matrix for EM
// doubtfull if omp is worth it
// matrix for nrm could be done much faster with a single call to elsym_i
// to check: dexter seems to use mean elsym
// expect long double to be sufficiently large and precise on all common archs

//[[Rcpp::export]]
arma::mat ittotmatC(const arma::vec& b, const arma::ivec& a, const arma::vec& c, arma::ivec& first, arma::ivec& last, const arma::ivec& ps)
{
	const int ms = accu(a(conv_to<uvec>::from(last)));
	const int nit = last.n_elem;
	const int npar = accu(last-first) + nit;
	const int nscores = ms+1;
	const vec logb = log(b);
	
	vec alogc(npar);
	for(int i=0;i<nit;i++)
		for (int j = first[i]; j <= last[i]; j++) 
			alogc[j] = a[j] * std::log(c[i]);

	  
	mat pi(npar, nscores, fill::zeros);	  
	
 
#pragma omp parallel
	{ 
		std::vector<long double> g(nscores), gw(nscores);
		vec eta(npar);
		int k;
		long double gs;
#pragma omp for	
		for (int s = 1; s <= ms; s++)
		{
			if(ps[s] == 1)
			{
				k = 0; // param
				eta = exp(logb + s * alogc);
				
				elsym_partial(eta,a, first.memptr(),last.memptr(),nit,g,gw,s,-1);
				gs = g[s];				
				
				for (int it = 0; it < nit; it++)
				{
					elsym_partial(eta,a, first.memptr(),last.memptr(),nit,g,gw,s-1,it);
					
					for (int j = first[it]; j <= last[it]; j++) 
					{
						if ( s-a[j] >= 0 && s-a[j] <= (ms - a[last[it]]) ) 
						{
							pi.at(k,s) = eta[j] * (g[s-a[j]]/gs);
						}
						k++;
					}
				}
			}
		}
	}
	return pi;
}



// matrix ga_s * gb_s / gab_s

// [[Rcpp::export]]
arma::mat sstable_nrmC(const arma::ivec& a, const arma::vec& b, 
						arma::ivec& firstA, arma::ivec& lastA, 
						arma::ivec& firstB, arma::ivec& lastB)
{
	const int ms_A = accu(a.elem(conv_to<uvec>::from(lastA)));
	const int ms_B = accu(a.elem(conv_to<uvec>::from(lastB)));
	const int nit_A = firstA.n_elem, nit_B = firstB.n_elem;
	
	ldvec gA(ms_A+1), gB(ms_B+1), g(ms_A + ms_B + 1);
	
	mat out(ms_A+1, ms_B+1,fill::zeros);
	
	elsym(b, a, firstA.memptr(), lastA.memptr(), nit_A, gA, g,-1);
	elsym(b, a, firstB.memptr(), lastB.memptr(), nit_B, gB, g,-1);
		
	g.zeros();
	for(int s=0; s<=ms_A; s++)
		for(int s2=0;s2<=ms_B;s2++)
			g[s+s2] += gA[s]*gB[s2];
			
	
	for(int s=0; s<=ms_A; s++)
		for(int s2=0;s2<=ms_B;s2++)
			out.at(s,s2) = gA[s]*gB[s2]/g[s+s2];
		

    return out;
}


// c must be per score, not per item
// [[Rcpp::export]]
arma::mat sstable_imC(const arma::ivec& a, const arma::vec& b, const arma::vec& c,					  
						arma::ivec& firstA, arma::ivec& lastA, 
						arma::ivec& firstB, arma::ivec& lastB)
{
	const int ms_A = accu(a.elem(conv_to<uvec>::from(lastA)));
	const int ms_B = accu(a.elem(conv_to<uvec>::from(lastB)));
	const int nit_A = firstA.n_elem, nit_B = firstB.n_elem;

	std::vector<long double> gA(ms_A+1), gB(ms_B+1), gw(std::max(ms_A,ms_B)+1);
	long double gs;
	int m1, start1;
	const vec logb = log(b), alogc = a % log(c);
	vec eta(b.n_elem);
	
	mat out(ms_A+1, ms_B+1,fill::zeros);	
	
	
	for(int s=1; s<=ms_A+ms_B; s++)
	{
		eta = exp(logb + s * alogc);	
		
		elsym_partial(eta,a, firstA.memptr(),lastA.memptr(),nit_A, gA,gw,s,-1);
		elsym_partial(eta,a, firstB.memptr(),lastB.memptr(),nit_B, gB,gw,s,-1);
		
		m1 = std::min(ms_A, s);
		start1 = std::max(0,s-ms_B);

		gs=0;
		for(int s1=start1; s1<=m1; s1++)
			gs += gA[s1] * gB[s-s1];
		
		for(int s1=start1; s1<=m1; s1++)
			out.at(s1,s-s1) = gA[s1]*gB[s-s1]/gs;

	}
	return out;
}






