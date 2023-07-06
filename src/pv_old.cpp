#include <RcppArmadillo.h>

using namespace arma;


double draw_theta(const vec& mu, const vec& sigma,  const double alpha)
{
	if(alpha>0)
	{
		int i = (R::runif(0,1) < alpha);
		return R::rnorm(mu[i],sigma[i]);
	}
	return R::rnorm(mu[0],sigma[0]);
}


// Produce PVs using recycling for normal or mixture normal

// scoretb is length maxscore + 1 and contains counts (possibly 0)
// theta's returned for persons sorted by score
// for simplicity, nPv is not taken into account. This can be done in R by multiplying scoretb by nPv and untying the result afterwards
// this has to be called separately per population (and per booklet)
// if alpha > 0, mu and sigma are assumed to have length 2 and denote a mixture normal
// otherwise length 1 and normal is assumed
// [[Rcpp::export]]
arma::vec PVrecycle(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last,
					const arma::vec& mu, const arma::vec& sigma, arma::ivec& scoretb, const arma::ivec& A, const double alpha = -1.0)
{
  
	int nP = accu(scoretb);
	const int nI = first.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const ivec cscoretb = cumsum(scoretb);
	  
	int x, k;
	double u, atheta;
	vec p(maxA+3, fill::zeros);  
	vec lookup(maxA+1);
	lookup[0] = 1.0;
	  
	vec theta(nP);
	  
	// observed maximum score
	int max_score = scoretb.n_elem - 1;  
	  
	for(;max_score>=0; max_score--)
		if(scoretb[max_score] > 0)
			break;

	int cntr = 0;
	while(nP > 0)
	{
		atheta = draw_theta(mu, sigma, alpha);
		for(int j=1;j<=maxA;j++) 
			lookup[j] = std::exp(j*atheta);
		x=0;
		for (int i=0;i<nI;i++)
		{
			p[0] = b[first[i]]; //  exp(0)==1
			k=1;
			for (int j=first[i]+1;j<=last[i];j++) // note the +1
			{
				p[k] = p[k-1] + b[j]*lookup[a[j]]; 
				k++;
			}
			u=p[k-1]*R::runif(0,1);
			k=0;
			while (u > p[k])
				k++;
			if (k > 0) 
				x += A[first[i]+k];
				
			if(x > max_score)
				break;
		}	
		if(scoretb[x] > 0)
		{
			theta[cscoretb[x] - scoretb[x]] = atheta;
			scoretb[x]--;
			nP--;

			if(x==max_score && scoretb[x] == 0)
				for(;max_score>=0; max_score--)
					if(scoretb[max_score] > 0)
						break;
		}
		if(++cntr > 50000)
		{
			Rcpp::checkUserInterrupt();
			cntr = 0;
		}		
	}
	return theta;
}


// [[Rcpp::export]]
void PV_slow(const arma::vec& b, const arma::ivec& a, const arma::ivec& bk_first, const arma::ivec& bk_last, 					
			const arma::ivec& bcni,
			const arma::ivec& booklet_id, const arma::ivec& booklet_score, const arma::vec& mu, const double sigma,
			arma::mat& pv_mat, const int pv_col_indx=0)
{
	const int np = pv_mat.n_rows;
	const int maxA = max(a);
	vec lookup(maxA+1);
	vec p(maxA+3, fill::zeros); 
	lookup[0] = 1.0;
	
	vec pv(pv_mat.colptr(pv_col_indx),np, false, true);
	
	double theta=0;
	int cntr=0;
	for(int prs=0; prs<np; prs++)
	{
		const int bk = booklet_id[prs];
		const int score = booklet_score[prs];
		
		int x=-1;
		while(x!=score)
		{		
			theta = R::rnorm(mu[prs], sigma);
			for(int j=1;j<=maxA;j++) 
				lookup[j] = std::exp(j*theta);
			
			x=0;
			for(int i=bcni[bk]; i<bcni[bk+1] && x <= score; i++)
			{
				p[0] = b[bk_first[i]]; 
				int k=1;
				for (int j=bk_first[i]+1;j<=bk_last[i];j++) 
				{
					p[k] = p[k-1] + b[j]*lookup[a[j]]; 
					k++;
				}
				double u=p[k-1]*R::runif(0,1);
				k=0;
				while (u > p[k])
					k++;
				x += a[bk_first[i]+k];
			}
			if(++cntr > 50000)
			{
				Rcpp::checkUserInterrupt();
				cntr = 0;
			}
		}
		pv[prs] = theta;

	}
}




// [[Rcpp::export]]
void PV_sve_old(const arma::vec& b, const arma::ivec& a, const arma::ivec& bk_first, const arma::ivec& bk_last, 					
			const arma::ivec& bcni,
			const arma::ivec& booklet_id, const arma::ivec& booklet_score, const arma::vec& mu, const double sigma,
			arma::mat& pv_mat, const int pv_col_indx=0, const int niter=1)
{
	const int np = pv_mat.n_rows;
	const int maxA = max(a);
	vec lookup(maxA+1);
	vec p(maxA+3, fill::zeros); 
	lookup[0] = 1.0;
	
	vec pv(pv_mat.colptr(pv_col_indx),np, false, true);
	
	for(int iter=0; iter<niter; iter++)
	{
		for(int prs=0; prs<np; prs++)
		{
			double theta = R::rnorm(mu[prs], sigma);
			for(int j=1;j<=maxA;j++) 
				lookup[j] = std::exp(j*theta);
			
			const int bk = booklet_id[prs];
			
			int x=0;
			for(int i=bcni[bk]; i<bcni[bk+1]; i++)
			{
				p[0] = b[bk_first[i]]; 
				int k=1;
				for (int j=bk_first[i]+1;j<=bk_last[i];j++) 
				{
					p[k] = p[k-1] + b[j]*lookup[a[j]]; 
					k++;
				}
				double u=p[k-1]*R::runif(0,1);
				k=0;
				while (u > p[k])
					k++;
				x += a[bk_first[i]+k];
			}
			
			double acc = std::exp((theta-pv[prs])*(booklet_score[prs]-x));
			if(R::runif(0,1)<acc)
				pv[prs] = theta;
		}
	}
}



