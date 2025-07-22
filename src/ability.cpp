#include <RcppArmadillo.h>
#include "ability.h"
#include "myomp.h"
#include "shared.h"
#include "elsym.h"



using namespace arma;

using Rcpp::Named;


// PA is working memory vector of length max_A + 1
template<bool LOGLIK>
void deriv_theta(const double theta, const vec& b, const ivec& a, int* first, int* last, const int nit, const int max_a, vec& PA, double& E, double& I, double &J)
{
	double M,M1,M2,M3,psum;
	
	E=0; I=0; J=0;
	
	for(int score=1; score<=max_a; score++)
		PA[score] = std::exp(score * theta);
	
	for(int i=0;i<nit;i++)
	{
		psum = 1;	// adaptation for omitting 0 cat
		for(int j = *(first+i); j<=*(last+i); j++)
			psum  += b[j] * PA[a[j]];

		M1=0; M2=0; M3=0;
		for(int j=*(first+i); j<=*(last+i); j++)
		{
			M = b[j] * PA[a[j]]/psum;	
			M1 += a[j] * M;
			M2 += SQR(a[j]) * M;
			M3 += CUB(a[j]) * M;			
		}
		if(LOGLIK)
		{
			E += std::log(psum);
		} else
		{
			E += M1;
		}			
		I += M2 - M1*M1;
		J += M3 - M1*(3.0*M2 - 2.0*M1*M1);
	} 
}

// MLE theta for a single score on an itemset (does not have to be an integer)
// [[Rcpp::export]]
arma::vec ML_theta_c(const double score, const arma::mat& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last)
{
	const int nit = first.n_elem, ndraws=b.n_cols, maxiter=200;
	const double acc = 1e-8;
	int max_a = 0;
	for(int i=0;i<nit;i++)
		max_a = std::max(max_a,a[last[i]]);
	
	vec wmem(max_a + 1), out_theta(ndraws);
	double E,I,J, theta=0;
	
	for(int draw=0;draw<ndraws;draw++)
	{
		deriv_theta<false>(theta, b.col(draw), a, first.memptr(), last.memptr(), nit, max_a, wmem, E, I, J);
		
		for(int iter=0; iter<maxiter; iter++)
		{
			E -= score;
			theta -= (2*E*I)/(2*SQR(I)-E*J);
			if(std::abs(E) < acc)
				break;
			deriv_theta<false>(theta, b.col(draw), a, first.memptr(), last.memptr(), nit, max_a, wmem, E, I, J);
		}
		out_theta[draw] = theta;	
	}
	return out_theta;
}

		


// [[Rcpp::export]]
Rcpp::List deriv_theta_c(const arma::vec& theta, const arma::mat& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last, const int n_cores=1)
{
	const int nit = first.n_elem, nt = theta.n_elem, ndraws=b.n_cols;
	int max_a = 0;
	for(int i=0;i<nit;i++)
		max_a = std::max(max_a,a[last[i]]);
	
	
	mat E(nt,ndraws), I(nt,ndraws), J(nt,ndraws);

#pragma omp parallel num_threads(n_cores)	
{
	vec wmem(max_a + 1);
#pragma omp for	
	for(int draw=0;draw<ndraws;draw++)
	{
		for(int t=0; t<nt; t++)
			deriv_theta<false>(theta[t], b.col(draw), a, first.memptr(), last.memptr(), nit, max_a, wmem, E.at(t,draw), I.at(t,draw), J.at(t,draw));
	}
}
	return Rcpp::List::create(Named("E")=E, Named("I")=I, Named("J")=J);
}

template<bool WLE>
double Escore_bk(const double theta, const vec& b, const ivec& a, int* first,  int* last, const int nit, const int max_a, vec& PA)
{
	double E,I,J;
	
	deriv_theta<false>(theta, b, a, first, last, nit, max_a, PA, E, I, J);
	
	if(WLE)
		E -= J/(2*I);
	
	return E;
}


Rcpp::List theta_output(mat& theta, mat& se, const ivec& bk_maxs, const ivec& bk_start, const int nbk, const bool add_inf=false)
{
	const int ndraws=theta.n_cols, nscores=theta.n_rows;
	
	ivec scores(nscores), booklet(nscores);	
	
	for(int bk=0;bk<nbk;bk++)
	{
		for(int s=0;s<=bk_maxs[bk];s++)
		{
			scores[s+bk_start[bk]] = s;
			booklet[s+bk_start[bk]] = bk+1; // R-indexed
		}
	}
	
	if(ndraws>1)
	{
		se.col(0) = sqrt(ndraws/(ndraws-1) * var(theta,0,1) + mean(square(se),1));		
		theta.col(0) = mean(theta,1);		
	} 
	
	if(add_inf)
	{
		for(int bk=0; bk<nbk; bk++)
		{
			theta.at(bk_start[bk]) = -1 * datum::inf;
			theta.at(bk_start[bk] + bk_maxs[bk],0) = datum::inf;
			se.at(bk_start[bk],0) = NA_REAL;
			se.at(bk_start[bk] + bk_maxs[bk], 0) = NA_REAL;
		}
	}
	return Rcpp::List::create(Named("booklet") = booklet, Named("booklet_score") = scores, Named("theta") = theta.col(0), Named("se") = se.col(0));
}


template<bool WLE>
Rcpp::List theta_wmle(const arma::mat& b, const arma::ivec& a, 
						arma::ivec& first, arma::ivec& last,
						const arma::ivec& bk_nit, const int n_cores=1)
{
	const int max_iter = 200;
	const double acc = 1e-8;
	const int ndraws = b.n_cols;
	
	const int nbk = bk_nit.n_elem;
	const ivec bk_cnit = ivec2_iter(bk_nit);
	
	ivec bk_maxs(nbk,fill::zeros), bk_max_A(nbk,fill::zeros);;
	
	for(int bk=0; bk<nbk; bk++)
	{
		for(int i=bk_cnit[bk]; i<bk_cnit[bk+1]; i++)
		{
			bk_maxs[bk] += a[last[i]];
			bk_max_A[bk] = std::max(bk_max_A[bk], a[last[i]]);
		}
	}
	const double max_A = max(bk_max_A);
	
	const ivec bk_cnscores = ivec2_iter(bk_maxs+1);
	const int nscores = accu(bk_maxs) + nbk;
	
	mat theta(nscores, ndraws, fill::zeros), se(nscores, ndraws, fill::zeros);
	
	// mle: for(int s=1; s<bk_maxs[bk]; s++)
	// wle: for(int s=0; s<=bk_maxs[bk]; s++)
	const int score_start = WLE ? 0 : 1;
	const int score_end = WLE ? 1 : 0;

	
#pragma omp parallel num_threads(n_cores)		
{
	double xl, rts, fl, f, dx;
	int *first_ptr, *last_ptr;
	
	bool conv_error;
	
	double E,I,J;
	vec wmem(max_A+1);
	
#pragma omp for
	for(int draw=0; draw<ndraws; draw++)
	{
		for(int bk=0; bk<nbk; bk++)
		{
			conv_error = false;
			xl = 0;
			rts = -1.3;
			
			first_ptr = first.memptr() + bk_cnit[bk];
			last_ptr = last.memptr() + bk_cnit[bk];
			
			fl = Escore_bk<WLE>(xl, b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem);
			f = Escore_bk<WLE>(rts, b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem);	

			// secant				
			for(int s=score_start; s<bk_maxs[bk]+score_end; s++)
			{
				for(int iter=0; iter<max_iter; iter++)
				{
					if(WLE && (xl > rts) != (fl > f))	conv_error = true;
				
					dx = (xl-rts) * (f-s)/(f-fl);
					xl = rts;
					fl = f;
					rts += std::copysign(std::min(std::abs(dx),0.5), dx); // steps larger than 0.5 on the theta scale are useless and can cause overflow in escore

					f = Escore_bk<WLE>(rts, b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem);
					
					if(std::abs(dx) < acc || conv_error)
						break;
				} 
				if(conv_error)
				{
					for(int s=score_start; s<bk_maxs[bk]+score_end; s++)
					{
						theta.at(s+bk_cnscores[bk],draw) = NA_REAL;
						se.at(s+bk_cnscores[bk],draw) = NA_REAL;
					}
					break;
				}
				theta.at(s+bk_cnscores[bk],draw) = rts;
				
				deriv_theta<false>(rts, b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem, E,I, J);
				
				if(WLE)
					se.at(s+bk_cnscores[bk],draw) = std::sqrt( (I+SQR(J/(2*I)))/SQR(I) );
				else
					se.at(s+bk_cnscores[bk],draw) = 1/std::sqrt(I);
				
				rts += 0.1; // give rts a nudge, otherwise (f-s)/(f-fl) can overflow since f-fl is often very small
				f = Escore_bk<WLE>(rts, b.col(draw), a, first_ptr, last_ptr, bk_nit[bk],  bk_max_A[bk], wmem);
			}
		}
	}
}

	return theta_output(theta, se, bk_maxs, bk_cnscores, nbk, !WLE);
}

// [[Rcpp::export]]
Rcpp::List theta_wmle_c(const arma::mat& b, const arma::ivec& a, 
						arma::ivec& first, arma::ivec& last,
						const arma::ivec& bk_nit, const bool WLE,
						const int n_cores=1)
{
	if(WLE)
		return theta_wmle<true>(b, a, first, last, bk_nit, n_cores);
	else
		return theta_wmle<false>(b, a, first, last, bk_nit, n_cores);
}



// eap jeffreys with a grid method

// [[Rcpp::export]]
Rcpp::List theta_jeap_c(const arma::vec& grid, const arma::mat& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last,  const arma::ivec& bk_nit, const int n_cores=1)
{
	// constants
	const int nt = grid.n_elem, nbk = bk_nit.n_elem, ndraws = b.n_cols;

	const ivec bk_cnit = ivec2_iter(bk_nit);	
	
	ivec bk_maxs(nbk,fill::zeros), bk_max_A(nbk,fill::zeros);
	
	for(int bk=0; bk<nbk; bk++)
	{
		for(int i=bk_cnit[bk]; i<bk_cnit[bk+1]; i++)
		{
			bk_maxs[bk] += a[last[i]];
			bk_max_A[bk] = std::max(bk_max_A[bk], a[last[i]]);
		}
	}
	
	const ivec bk_start = ivec2_iter(bk_maxs+1);
	
	const int nscores = accu(bk_maxs) + nbk;
	
	//output vars
	mat theta(nscores,ndraws,fill::zeros), se(nscores,ndraws,fill::zeros);	
	
#pragma omp parallel num_threads(n_cores)		
{	
	// working variables
	int *first_ptr, *last_ptr;
	double E,I,J;
	long double w, sum_w, sum_theta;
	vec wmem(max(bk_max_A)+1);
	vec prior(nt,fill::zeros), ll(nt,fill::zeros);
	
#pragma omp for
	for(int draw=0; draw<ndraws; draw++)
	{	
		for(int bk=0; bk<nbk; bk++)
		{
			first_ptr = first.memptr() + bk_cnit[bk];
			last_ptr = last.memptr() + bk_cnit[bk];
		
			for(int t=0;t<nt;t++)
				deriv_theta<true>(grid[t], b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem, ll[t], prior[t], J);			

			prior = sqrt(prior);
			
			for(int s=0; s<=bk_maxs[bk]; s++)
			{
				sum_w=0, sum_theta=0;
				for(int t=0;t<nt;t++)
				{
					w = prior[t] * std::exp(s * grid[t] - ll[t]);
					sum_w += w;
					sum_theta += w * grid[t];
				}
				theta.at(bk_start[bk]+s, draw) = sum_theta/sum_w;
				deriv_theta<false>(theta.at(bk_start[bk]+s, draw), b.col(draw), a, first_ptr, last_ptr, bk_nit[bk], bk_max_A[bk], wmem, E, I, J);	
				se.at(bk_start[bk]+s, draw) = std::sqrt( (I+SQR(J/(2*I))) /SQR(I));
			}
		}
	}
}
	return theta_output(theta, se, bk_maxs, bk_start, nbk, false);

}

// eap normal with a grid

double logsumexp(const vec& v)
{
	const double m = max(v);
	return m + std::log(accu(exp(v-m)));
}


// [[Rcpp::export]]
Rcpp::List theta_eap_c(const arma::vec& grid, const arma::vec& weights, const arma::mat& b, const arma::ivec& a, arma::ivec& first, arma::ivec& last,  const arma::ivec& bk_nit, const int n_cores=1)
{
	// constants
	const int nt = grid.n_elem, nbk = bk_nit.n_elem, ndraws = b.n_cols;

	const ivec bk_cnit = ivec2_iter(bk_nit);	
	
	ivec bk_maxs(nbk,fill::zeros);
	
	for(int bk=0; bk<nbk; bk++)
		for(int i=bk_cnit[bk]; i<bk_cnit[bk+1]; i++)
			bk_maxs[bk] += a[last[i]];

	const ivec bk_start = ivec2_iter(bk_maxs+1);
	
	const int nscores = accu(bk_maxs) + nbk, max_score = max(bk_maxs);
	
	//output vars
	mat theta(nscores,ndraws,fill::zeros), se(nscores,ndraws,fill::zeros);		
	
	
#pragma omp parallel num_threads(n_cores)		
{	
	// working variables
	int *first_ptr, *last_ptr;
	long double sumw, m;
	vec g(max_score+1), gw(max_score+1), w(nt);
	mat ps(max_score+1, nt);
	
	
#pragma omp for
	for(int draw=0; draw<ndraws; draw++)
	{		
		for(int bk=0; bk<nbk; bk++)
		{
			first_ptr = first.memptr() + bk_cnit(bk);
			last_ptr = last.memptr() + bk_cnit(bk);
			elsym(b.col(draw), a, first_ptr, last_ptr, bk_nit(bk), g, gw);
			g = log(g);

			for(int t=0; t<nt; t++)
			{
				for(int s=0; s<=bk_maxs[bk]; s++)
					ps.at(s,t) = g[s] + s * grid[t];
				ps.col(t).head(bk_maxs[bk]+1) -= logsumexp(ps.col(t).head(bk_maxs[bk]+1));
			}
			ps = exp(ps); // probability score s 

			for(int s=0; s<=bk_maxs[bk]; s++)
			{
				sumw=0;
				for(int t=0; t<nt; t++)
				{
					w[t] = ps.at(s,t) * weights[t];
					sumw += w[t]; 
				}
				w /= sumw;
				
				m = accu(grid % w);
				theta.at(bk_start[bk]+s, draw) = m;
				se.at(bk_start[bk]+s, draw) = std::sqrt(accu(square(grid-m) % w));	
			}			
		}
	}
}	
	return theta_output(theta, se, bk_maxs, bk_start, nbk, false);
}	
