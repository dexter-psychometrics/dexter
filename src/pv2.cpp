#include <RcppArmadillo.h>

using namespace arma;

using Rcpp::Named;

// two groups, one booklet

// [[Rcpp::export]]
arma::imat score_tab_double(const arma::ivec& scores, const arma::ivec& group, const int max_score)
{
	const int n = scores.n_elem;
	imat out(max_score+1, 2, fill::zeros);
	
	for(int i=0; i<n; i++)
		out.at(scores[i], group[i])++;
 
	return out;
}

// sort of the reverse operation of the previous to get the theta's in the proper order

// [[Rcpp::export]]
arma::vec arrange_pv(const arma::vec& theta, const arma::ivec& scoretab_np, const arma::ivec& group)
{
	vec out(theta.n_elem, fill::zeros); // fill zeros to test, is not necessary
	
	const int ntab = scoretab_np.n_elem;
	
	ivec cntr(2, fill::zeros);	
	
	int pstart = 0; 
	
	for(int tab=0; tab < ntab; tab+=2)
	{
		cntr[0] = pstart;
		cntr[1] = pstart + scoretab_np[tab];
		
		const int pend = pstart + scoretab_np[tab] + scoretab_np[tab+1];
		
		for(int p=pstart; p<pend; p++)
		{
			out[p] = theta(cntr[group[p]]++);
		}
		pstart = pend;
	}
	return out;
}

// omitted for now: user interrupt check
// omitted alpha because it was never used
// wil make templated for alpha possibly, or a class for draw when we no the parallel, mu and sigma can be supplied as 1 or 2 column matrices

// [[Rcpp::export]]
void pv_draw(const arma::vec& b, const arma::ivec& a, const arma::ivec& A, const arma::ivec& first, const arma::ivec& last, 
			 const arma::ivec& bk_cnit, const arma::ivec& bk_max_a,
		     const arma::ivec& const_scoretab, const arma::ivec& cscoretab, const arma::ivec& scoretab_bk, const arma::ivec& scoretab_pop, 
			 const arma::ivec& scoretab_nscores, const arma::ivec& scoretab_cnscores, const arma::ivec& scoretab_np,
			 const arma::vec& mu, const arma::vec& sigma,
			 arma::mat& theta, const int col_indx)
{
	const int ntab = scoretab_bk.n_elem;
	
	ivec scoretab = const_scoretab; // will get eaten

	// parallel init
	vec expat(bk_max_a.max()+1, fill::zeros), p(bk_max_a.max()+1, fill::zeros);
	
	int x,y,k;
	double u, atheta;
	
	// parallel for
	for(int tab=0; tab < ntab; tab++)
	{
		int np = scoretab_np[tab];
		int max_score = scoretab_nscores[tab] - 1;			
		
		const int scoretab_start = scoretab_cnscores[tab];
		
		const int bk = scoretab_bk[tab];
		const int item_start = bk_cnit[bk];
		const int item_end = bk_cnit[bk+1];
		const int max_a = bk_max_a[bk];
		
		const double pmu = mu[scoretab_pop[tab]];
		const double psigma = sigma[scoretab_pop[tab]];		
		
		expat[0] = 1.0;
		
		for(;max_score>=0; max_score--)
			if(scoretab[scoretab_start + max_score] > 0)
				break;
		
		while(np>0)
		{
			atheta = R::rnorm(pmu, psigma);
			
			for(int j=1; j<=max_a; j++) 
				expat[j] = std::exp(j*atheta);
			
			x = 0;
			for (int i=item_start; i<item_end; i++)
			{
				p[0] = b[first[i]]; 
				k=1;
				for (int j=first[i]+1;j<=last[i];j++) // note the +1
				{
					p[k] = p[k-1] + b[j] * expat[a[j]]; 
					k++;
				}
				u = p[k-1]*R::runif(0,1);
				k=0;
				while (u > p[k])
					k++;
				if (k > 0) // if seems unnecessary since zero cat should be in a and A?
					x += A[first[i]+k];
					
				if(x > max_score)
					break;
			}	
			y = scoretab_start + x;
			
			if(scoretab[y] > 0)
			{
				theta.at(cscoretab[y] - scoretab[y], col_indx) = atheta;
				scoretab[y]--;
				np--;

				if(x == max_score && scoretab[y] == 0)
					for(;max_score>=0; max_score--)
						if(scoretab[scoretab_start + max_score] > 0)
							break;
			}		
		}	
	}
	//printf("scoretab: min: %i, max: %i, 0: %i\n",scoretab.min(),scoretab.max(),accu(scoretab));
	//fflush(stdout);
}



