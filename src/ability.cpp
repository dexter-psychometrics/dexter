#include <RcppArmadillo.h>
#include "ability.h"

using namespace arma;




/*****
* vectorized version of Escore
*
* optimized with max_a
*****/

// [[Rcpp::export]]
arma::vec Escore_C(const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
				 const arma::ivec& first, const arma::ivec& last)
{
	const int nI = first.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int m = theta.n_elem;  

	double denom, num;
	vec lookup(maxA+1);
	
	vec score(m, fill::zeros);

	lookup[0] = 1.0;
	  
	for(int p=0; p<m; p++)
	{
		for(int i=1; i <= maxA; i++)
		{
			lookup[i] = exp(i*theta[p]);
		}
	  
		for (int i=0; i<nI; i++)
		{
			num=0.0;
			denom=1.0;
			for (int j=first[i]+1;j<=last[i];j++) // note +1
			{
			  num    +=a[j]*b[j]*lookup[a[j]];
			  denom  +=     b[j]*lookup[a[j]];
			}
			score[p]+=num/denom;
		}
	}
	return score;
}




// non vectorized for use in  theta_mle
double Escore_single(double theta, const vec& b, const ivec& a, const ivec& first,  const ivec& last, const int n, const int max_a)
{

  double denom, num;
  vec lookup(max_a+1);
  double score = 0;

  lookup[0] = 1.0;
  

  for(int i=1; i <= max_a; i++)
  {
    lookup[i] = exp(i*theta);
  }
  
  for (int i=0; i<n; i++)
  {
    num=0.0;
    denom=1.0;
    for (int j=first[i]+1;j<=last[i];j++) // note +1
    {
      num    +=a[j]*b[j]*lookup[a[j]];
	  denom  +=     b[j]*lookup[a[j]];
    }
    score+=num/denom;
  }
  return score;
}




// Secant method used to find MLE of theta 
// [[Rcpp::export]]
arma::vec theta_mle_sec(const arma::vec& b, const arma::ivec& a, 
						const arma::ivec& first, const arma::ivec& last)
{
	const int n = first.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int max_score = accu(a(conv_to<uvec>::from(last)));
	
	vec theta(max_score-1);
	
	double xl = 0, rts = -1.3;
	double fl = Escore_single(xl, b, a, first, last, n, maxA),
		   f = Escore_single(rts, b, a, first, last, n, maxA);
	
	double dx;
	
	const int max_iter = 200;
	const double acc = 1e-8; // this might be a little too small but it seems to work

	// secant	
	for(int s=1; s<max_score; s++)
	{
		for(int iter=0; iter<max_iter; iter++)
		{
			dx = (xl-rts) * (f-s)/(f-fl);
			xl = rts;
			fl = f;
			rts += std::copysign(std::min(std::abs(dx),0.5), dx); // steps larger than 0.5 on the theta scale are useless and can cause overflow in escore
			f = Escore_single(rts, b, a, first, last, n, maxA);
			
			if(std::abs(dx) < acc)
				break;
		} 
		theta[s-1] = rts;
		rts += 0.1; // give rts a nudge, otherwise (f-s)/(f-fl) can overflow since f-fl is often very small
		f = Escore_single(rts, b, a, first, last, n, maxA);
	}
	return theta;
}


// [[Rcpp::export]]
double escore_wle(const double theta, const arma::vec& b, const arma::ivec& a, const arma::ivec& first,  const arma::ivec& last, const int nI, const int max_a)
{
	const int max_ncat = max(last - first) + 1;
	
	std::vector<long double> Fij(max_ncat);
	
	long double I=0,J=0;


	for(int i=0;i<nI;i++)
	{
		long double colsm = 0;	
		for(int j = first[i],k=0; j<=last[i]; j++, k++)
		{
			Fij[k] = b[j] * exp(a[j] * theta);
			colsm += Fij[k];
		}

		long double M1=0, M2=0, M3=0;
		for(int j=first[i], k=0; j<=last[i];j++, k++)
		{
			long double M = Fij[k]/colsm;	
			M1 += a[j] * M;
			M2 += (a[j] * a[j]) * M;
			M3 += (a[j] * a[j] * a[j]) * M;			
		}
			
		I += M2 - M1*M1;
		J += M3 - M1*(3.0*M2 - 2.0*M1*M1);
	}
	
	return Escore_single(theta, b, a, first, last, nI, max_a) - (J/(2*I));
}



// Secant method used to find W-MLE of theta
// gives an error for non-monotonicity 
// [[Rcpp::export]]
arma::vec theta_wle_sec(const arma::vec& b, const arma::ivec& a, 
						const arma::ivec& first, const arma::ivec& last)
{
	const int n = first.n_elem;
	const int maxA = max(a(conv_to<uvec>::from(last))); 
	const int max_score = accu(a(conv_to<uvec>::from(last)));
	
	vec theta(max_score+1);
	
	double xl = 0, rts = -1.3;
	double fl = escore_wle(xl, b, a, first, last, n, maxA),
		   f = escore_wle(rts, b, a, first, last, n, maxA);
	
	double dx;
	
	const int max_iter = 200;
	const double acc = 1e-8; // this might be a little too small but it seems to work

	// secant	
	for(int s=0; s<=max_score; s++)
	{
		for(int iter=0; iter<max_iter; iter++)
		{
			// protection against non-monotonicity
			if((xl > rts) != (fl > f))
			{
				//printf("score: %i\nxl: %f > rts: %f\nfl: %f > f %f\n", s, xl, rts, fl, f);
				//fflush(stdout);
				Rcpp::stop("Warm WLE estimates do not converge");
			}
			
			dx = (xl-rts) * (f-s)/(f-fl);
			xl = rts;
			fl = f;
			//rts += dx;
			rts += std::copysign(std::min(std::abs(dx),0.5), dx); 
			f = escore_wle(rts, b, a, first, last, n, maxA);						
			if(std::abs(dx) < acc)
				break;
			
		} 
		theta[s] = rts;
		rts += 0.1; // give rts a nudge, otherwise (f-s)/(f-fl) can overflow since f-fl is often very small
		f = escore_wle(rts, b, a, first, last, n, maxA);
	}
	return theta;
}

//only reason for this to be in c is extended precision type
// [[Rcpp::export]]
Rcpp::List theta_EAP_GH_c(const arma::mat& p_score, const arma::vec& theta, const arma::vec& weights)
{
	const int ns = p_score.n_cols, nt=p_score.n_rows;
	std::vector<long double> w(nt);
	vec eap(ns), se(ns);
	for(int s=0; s<ns; s++)
	{
		long double sumw=0;
		for(int t=0; t<nt; t++)
		{
			w[t] = p_score.at(t,s) * weights[t];
			sumw += w[t]; 
		}
		long double m=0, var=0;
		for(int t=0; t<nt; t++)	m += theta[t] * w[t]/sumw;
		for(int t=0; t<nt; t++)	var += (theta[t] - m) * (theta[t] - m) * w[t]/sumw;
		eap[s] = m;
		se[s] = std::sqrt(var);
	}
	return Rcpp::List::create(Rcpp::Named("theta") = eap, Rcpp::Named("se") = se);
}


/*
* output: I, J, logFi assumed to be initialized to zero by caller
*
* probably still not an optimal implementation but close
*/
// [[Rcpp::export]]
void IJ_c(const arma::vec& theta, const arma::vec& b, const arma::ivec& a, 
		  const arma::ivec& first, const arma::ivec& last, 
		  arma::mat& I, arma::mat& J,  arma::vec& logFi)
{
	const int nI = first.n_elem;
	const int nT = theta.n_elem;	
	const int max_ncat = max(last - first) + 1;
	
	mat Fij(nT, max_ncat);
	double M,M1,M2,M3,colsm;

	const ivec a2 = a % a;
	const ivec a3 = a2 % a;	
	
	for(int i=0;i<nI;i++)
	{
		for(int k=0;k<nT;k++)
		{
			for(int j = first[i]; j<=last[i]; j++)
			{
				Fij.at(k, j-first[i]) = b[j] * exp(a[j] * theta[k]);
			}
		}

		for(int j=0;j<nT;j++)
		{
			colsm = 0.0;	
			for(int q=0;q<=last[i] - first[i];q++)
			{
				colsm  += Fij.at(j,q);
			}
			logFi[j] += log(colsm);
			
			M1=0;M2=0;M3=0;
			for(int k=first[i]; k<=last[i];k++)
			{
				M = Fij.at(j, k-first[i])/colsm;	
				M1 += a[k] * M;
				M2 += a2[k] * M;
				M3 += a3[k] * M;			
			}
			
			I.at(j,i) = M2 - M1*M1;
			J.at(j,i) = M3 - M1*(3.0*M2 - 2.0*M1*M1);
		}

	} 
}

