#include <RcppArmadillo.h>
#include <xoshiro.h>
#include <omp.h>
#include "progress.h"
#include "priors.h"


using namespace arma;

#define SEED std::round(R::runif(0,1) * 2147483647)

void ElSym(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS,  
			std::vector<long double>& g) 
{
  int Msc = 0, idx = 0; // start from column one
  
  std::vector<long double> gg(2 * mS + 2, 0);
  std::fill(g.begin(),g.end(),0);
  gg[0] = 1;

  for (int i=0; i<nI; i++)
  {
	if (i != item1 && i != item2)
    {
      for (int s=Msc; s>=0; s--)
		gg[2*s+(1-idx)] = 0;
      
	  for (int s=Msc; s>=0; s--)
        for (int j=first[i]; j<=last[i]; j++)
          if (b[j]>0)
			gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]; 

      Msc += a[last[i]];
      idx = 1-idx; // swap columns
    }
  }
  
  for (int s=0; s<=Msc; s++)
	g[s] = gg[2*s + idx];

}


//[[Rcpp::export]]
void ElSym_C(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, arma::vec& g)
{
	std::vector<long double> gl(g.n_elem);
	ElSym(b,a,first,last,item1, item2,nI,mS,gl);
	g.zeros();
	for(int s=0;s<=mS;s++)
		g[s] = gl[s];
} 

// compute the gamma functions g1(-item1) and g2(-item1, -item2)
void ElSym12(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, 
				std::vector<long double>& g1, std::vector<long double>& g2) 
{
	int Msc = 0, idx = 0; // start from column one
	
	  
	std::vector<long double> gg(2 * mS + 2,0);
	std::fill(g1.begin(),g1.end(),0);
	std::fill(g2.begin(),g2.end(),0);

	gg[0] = 1;

	for (int i=0; i<nI; i++)
	{
		if (i != item1 && i != item2)
		{
			for (int s=Msc; s>=0; s--)
				gg[2*s+(1-idx)] = 0;
		  
			for (int s=Msc; s>=0; s--)
				for (int j=first[i]; j<=last[i]; j++)
					if (b[j]>0)
						gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]; 

			Msc += a[last[i]];
			idx = 1-idx; // swap columns
		}
	}
	// store g2  
	for (int s=0; s<=Msc; s++)
		g2[s] = gg[2*s + idx];
	
	// 1 loop remaining to add item2 to g1
	for (int s=Msc; s>=0; s--)
		gg[2*s+(1-idx)] = 0;
	
	for (int s=Msc; s>=0; s--)
		for (int j=first[item2]; j<=last[item2]; j++)
			if (b[j]>0)
				gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]; 

	Msc += a[last[item2]];
	idx = 1-idx;
	
	// store g1  
	for (int s=0; s<=Msc; s++)
		g1[s] = gg[2*s + idx];
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
				for(int j=first[i]+1; j<=last[i]; j++)
					ps(s+a[j]) = 1;
		m += a[last[i]];
	}	
	return ps;
}





//[[Rcpp::export]]
arma::mat ittotmat0_C(const arma::vec& b, const arma::ivec& a, const arma::vec& c, const arma::ivec& first, const arma::ivec& last, const arma::ivec& ps)
{
	const int ms = accu(a(conv_to<uvec>::from(last)));
	const int nI = last.n_elem;
	const int npar = accu(last-first) + nI;
	const int nscores = ms+1;
	const int i1=-1, i2=-1;  
	const vec logb = log(b);
	const vec alogc = a % log(c);  
	  
	mat pi(npar, nscores, fill::zeros);	  
 
#pragma omp parallel
	{ 
		std::vector<long double> g(nscores), gi(nscores);
		vec eta(npar);
#pragma omp for	
		for (int s = 0; s <= ms; s++)
		{
			if(ps[s] == 1)
			{
				int k = 0; 
				eta = exp(logb + s * alogc);
				
				ElSym(eta, a, first, last, i1, i2, nI, ms, g);
				
				for (int it = 0; it < nI; it++)
				{
					ElSym(eta, a, first, last,  it, i2, nI, ms, gi);
					for (int j = first[it]; j <= last[it]; j++) 
					{
						int indx = s-a[j];
						if ( indx >= 0 && indx < (nscores - a[last[it]]) ) 
						{
							//pi.at(k,s) = exp(log(eta[j]) + log(gi[indx]) - log(g[s]));
							pi.at(k,s) = eta[j] * (gi[indx]/g[s]);
						}
						k++;
					}
				}
			}
		}
	}
	return pi;
}





/*
computes Hessian and E.step/Expect as collateral, saving some calls to Elsym
b,a,first, last are presumed to refer to a single booklet

*/

void E_Hess(const vec& b, const ivec& a, const ivec& first, const ivec& last, int* ptr_scoretab, const int nI, const int n_score, vec& expect, mat& Hess)
{

	const int i1=-1, i2=-1;
	const int nPar = b.n_elem;
	const int mS = n_score - 1;
	double tmp;
	vec cc(mS+1, fill::zeros);
	
	const ivec scoretab(ptr_scoretab, n_score, false, true);	  
	  
	std::vector<long double> g(mS+1), gi(mS+1), gk(mS+1), gik(mS+1);
	
	ElSym(b,a,first,last,i1,i2,nI,mS,g);

	for (int item=0;item<nI;item++)
	{
		ElSym(b,a,first,last,item,i2,nI,mS,gi);
		for (int j=(first[item]+1);j<=last[item];j++)
		{
		  for (int s=a[j];s<=mS;s++)
		  {		    
			if (g[s]>0)
			{
			  tmp = b[j] * ( gi[s-a[j]] / g[s] );
			  cc[s] = scoretab[s] * tmp;		
			  expect[j] += cc[s]; // this line replaces function Expect
			  Hess.at(j,j) += cc[s] * (1-tmp); // within current item
			  
			}
		  }
		  for (int k=(j+1);k<=last[item];k++)
		  {
			for (int s=a[k];s<=mS;s++)
			{
			  if (g[s]>0)
			  {
				Hess.at(k,j) -= cc[s] * b[k] * (gi[s-a[k]]/g[s]); //within current item 
			  }
			}
		  }

		  for (int k=(item+1);k<nI;k++) //deze loopt over alle volgende items
		  {
			//ElSym(b,a,first,last,k,i2,nI,mS,gk);
			//ElSym(b,a,first,last,k,item,nI,mS,gik);
			ElSym12(b,a,first,last,k,item,nI,mS,gk,gik);

			for (int l=(first[k]+1);l<=last[k];l++)
			{
			  for (int s=0;s<=mS;s++)
			  {
				if (g[s]>0)
				{
				  if (s >= (a[j]+a[l]))
				  { 
					Hess.at(l,j) += scoretab[s] * b[j] * b[l] * (gik[s-a[j]-a[l]] / g[s]); //j is current item, l is andere(opvolgende) items
				  }
				  if ((s>=a[j])&&(s>=a[l])) 
				  { 
					Hess.at(l,j) -= cc[s] * b[l] * (gk[s-a[l]] / g[s]);
				  }               
				}
			  }
			}
		  }
		}
	}

	// sum Hess with transpose excepting diagonal;
	for(int i=0; i<nPar; i++)
		for(int j=i+1; j<nPar; j++)
		{
			tmp = Hess.at(i,j);
			Hess.at(i,j) += Hess.at(j,i);
			Hess.at(j,i) += tmp;			
		}
	
}


void Expect(const vec& b, const ivec& a, int* ptr_first, int* ptr_last, int* ptr_scoretab, int nI, int n_score, vec& expect)
{
  const int i1=-1, i2=-1;
  const int mS = n_score - 1;
  const ivec first(ptr_first, nI, false, true), last(ptr_last, nI, false, true);
  const ivec scoretab(ptr_scoretab, n_score, false, true);	 

  
  std::vector<long double> g(n_score, 0), gi(n_score, 0);
  
  ElSym(b,a,first,last,i1,i2,nI,mS,g);
  
  for (int item=0; item<nI; item++)
  {    
	ElSym(b,a,first,last,item,i2,nI,mS,gi);
    for (int j=first[item];j<=last[item];j++)
      for (int s=a[j];s<=mS;s++)
        if (g[s]>0) 
			expect[j] += scoretab[s] * b[j] * ( gi[s-a[j]] / g[s] );
  }
}



/* is parallellizable but not much to gain 
*
* done in this matter to have call convention similar to NR_booklets
* which saves code in anon
*/
// [[Rcpp::export]]
arma::vec E_booklets(const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit)
{
	const int nb = nit.n_elem;
	const int nPar = b.n_elem;
	vec E(nPar, fill::zeros);
	vec Et(nPar);
	
	// cumulatives for bookkeeping
	ivec cnit(nb), cn_score(nb); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nb-1) = cumsum(nit.head(nb-1));
	cn_score.tail(nb-1) = cumsum(n_score.head(nb-1));
	
	for(int bl=0; bl < nb; bl++)
	{
		Et.zeros();
		Expect(b, a, first.memptr() + cnit(bl), last.memptr() + cnit(bl), scoretab.memptr() + cn_score(bl), nit(bl), n_score(bl), Et);
		E += Et;
	}
	return E;
}



/*
first and last need to be zero indexed!
do E_hess for all booklets with omp optimization

made it subroutine to easily return both esufi and hess

everything except b and a is concatenated per booklet
*/

// to do: use reductions in omp, code can be done much simpler

// [[Rcpp::export]]
void NR_booklets(const arma::vec& b, const arma::ivec& a, 
				 arma::ivec& first, arma::ivec& last, arma::ivec& scoretab, 
				 const arma::ivec& n_score, const arma::ivec& nit, const int max_par_bk,
				 arma::vec& EsufI, arma::mat& H)
{
	const int nb = nit.n_elem;
	//const int nPar = EsufI.n_elem;
	const int max_nit_bk = max(nit);
	
	EsufI.zeros();
	H.zeros();
	
	// if just one booklet, avoid unnecessary overhead in parallelisation
	if(nb==1)
	{
		E_Hess(b, a, first, last, scoretab.memptr(), nit(0), n_score(0), EsufI, H);
		return;
	}
	
	// cumulatives for bookkeeping
	ivec cnit(nb), cn_score(nb); 
	
	cnit[0] = 0;
	cn_score[0] = 0;
	cnit.tail(nb-1) = cumsum(nit.head(nb-1));
	cn_score.tail(nb-1) = cumsum(n_score.head(nb-1));
	

#pragma omp parallel
	{ 
		mat Ht(max_par_bk, max_par_bk);
		vec Et(max_par_bk), bkb(max_par_bk);
		ivec bkfirst(max_nit_bk), bklast(max_nit_bk), bka(max_par_bk);
		
#pragma omp for			
		for(int bl=0; bl < nb; bl++)
		{
			Ht.zeros();
			Et.zeros();
			// easier to make cnit 1 longer
			// bookkeeping
			int itm=0, par=0;
			for(int i=cnit[bl]; i< (cnit[bl] + nit[bl]); i++)
			{
				bkfirst[itm] = par;				
				for(int j=first[i];j<=last[i];j++)
				{
					bka[par] = a[j];
					bkb[par++] = b[j];
				}
				bklast[itm++] = par-1;
			}
	
			E_Hess(bkb, bka, bkfirst, bklast, scoretab.memptr() + cn_score(bl), nit(bl), n_score(bl), Et, Ht);
		
#pragma omp critical
			{
				int p=0, q;
				for(int i=cnit[bl]; i< (cnit[bl] + nit[bl]); i++)
				{
					for(int j=first[i];j<=last[i];j++)
					{
						EsufI[j] += Et[p];
						q=0;
						for(int i2=cnit[bl]; i2<(cnit[bl] + nit[bl]); i2++)	
							for(int j2=first[i2]; j2<=last[i2]; j2++)
								H.at(j,j2) += Ht.at(p, q++);
						p++;
					}
				}
			}	
		}
	}
}





// Bayesian calibration

void elsym_helper(const vec& b, const ivec& a, int* ptr_first, int* ptr_last, const int item1, const int nI, const int mS, std::vector<long double>& g)
{
	const ivec first(ptr_first, nI, false, true), last(ptr_last, nI, false, true);
	
	ElSym(b, a, first, last, item1, -1, nI, mS, g);
}



// [[Rcpp::export]]
Rcpp::List calibrate_Bayes_C(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, 
				 const arma::ivec& ib, const arma::ivec& bi, const arma::ivec& nbi, 
				 const arma::ivec& nib,
				 arma::ivec& bfirst, arma::ivec& blast,
				 const arma::ivec& bmax, const arma::ivec& m,
				 const arma::ivec& sufI, const arma::ivec& bkscoretab,
				 const arma::vec& b_in,  const arma::vec& fixed_b, 
				 const int from, const int step, const int ndraws,
				 const arma::ivec progress_init,
				 const double prior_eta=0.5, const double prior_rho=0.5)
{

	const bool free_calibration = all(fixed_b != fixed_b); // NA != NA
	const int nIter = ndraws * step + from;
	const int max_cat = max(last-first)+1;
	
	progress pb(ndraws, progress_init);
	
	// cumulatives for bookkeeping
	const int nI = nbi.n_elem;
	ivec cnbi(nI+1);
	cnbi(0) = 0;
	cnbi.tail(nI) = cumsum(nbi);
	
	const int nB = nib.n_elem;
	ivec cnib(nB+1);
	cnib(0) = 0;
	cnib.tail(nB) = cumsum(nib);
	
	ivec cbmax(nB+1);
	cbmax[0] = 0;	
	for(int k=0; k<nB; k++)
		cbmax[k+1] = cbmax[k] + bmax[k]+1;
	

	const int max_bscore = max(bmax);

	// working variables
	vec y(max_cat, fill::zeros), z(nB, fill::zeros);
	vec bklambda(bkscoretab.n_elem, fill::ones);	// to do: fill 0 lijkt beter, onmogelijke scores?
	mat bklambdax(ndraws, bkscoretab.n_elem);
	
	std::vector<long double> g(max_bscore+1);
	vec pi_k(max_bscore+1, fill::zeros);
	//vec fpwr(max_bscore+1);
	//fpwr[0] = 1;
	
	// copy b
	vec b(b_in.memptr(), b_in.n_elem);
	mat bx(ndraws, b.n_elem);
	
	//Gibbs
	int k,row_index=0;
	long double sm;
	
	for (int iter=0; iter<nIter; iter++)
	{
		for (int k=0; k<nB; k++)
		{
			elsym_helper(b,a, bfirst.memptr() + cnib[k], blast.memptr() + cnib[k], -1, nib[k], bmax[k], g);
			
			sm = 0;
			for (int s=0; s<=bmax[k];s++)
			{
			  if (g[s]>0){
			    pi_k[s] = R::rgamma(bkscoretab[cbmax[k]+s]+0.1, 1.0);
			    sm += pi_k[s];
			  }
			}
			
			for (int s=0; s<=bmax[k];s++)
			{
			  if(g[s]>0)
			  {
			    pi_k[s] = pi_k[s]/sm;
			    bklambda[cbmax[k]+s] = (pi_k[s]*m[k])/g[s];
			  }
			}
			
			z[k] = R::rgamma(m[k], 1.0/m[k]);
			
		}
		// bayes_items
		for (int i=0; i<nI;i++)
		{
			y.zeros();
			for (int ib_indx=cnbi[i]; ib_indx<cnbi[i+1]; ib_indx++)
			{
				k = ib[ib_indx];
				elsym_helper(b, a, bfirst.memptr() + cnib[k], blast.memptr() + cnib[k], bi[ib_indx], nib[k], bmax[k], g);
					
				for (int s=0; s<=bmax[k]-a[last[i]]; s++)
					for (int j=first[i]+1, c=0; j<=last[i]; j++,c++)
						y[c] += z[k] * g[s] * bklambda[cbmax[k]+s+a[j]];			
			}

			for (int j=first[i]+1, c=0; j<=last[i]; j++,c++)
				b[j] = R::rgamma(sufI[j] + prior_eta, 1/(y[c] + prior_rho));
		}
		

		if (!free_calibration)
		{
			b.elem(find_finite(fixed_b)) = fixed_b.elem(find_finite(fixed_b));
		}

		b.replace(datum::nan, 1); 
		if(iter >= from && (iter-from) % step == 0)
		{
			Rcpp::checkUserInterrupt();
			bklambdax.row(row_index) = bklambda.t();
			bx.row(row_index++) = b.t();
			pb.tick();
		}
		
	}
	return Rcpp::List::create(Rcpp::Named("b") = bx, Rcpp::Named("lambda") = bklambdax);
}


// [[Rcpp::export]]
Rcpp::List calibrate_Bayes_chains(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, 
				 const arma::ivec& ib, const arma::ivec& bi, const arma::ivec& nbi, 
				 const arma::ivec& nib,
				 arma::ivec& bfirst, arma::ivec& blast,
				 const arma::ivec& bmax, const arma::ivec& m,
				 const arma::ivec& sufI, const arma::ivec& bkscoretab,
				 const arma::mat& b_start,  const arma::vec& fixed_b, 
				 const int warmup, const int step, const int ndraws,
				 const arma::ivec progress_init, const int max_cores,
				 const double prior_eta=0.5, const double prior_rho=0.5)
{
	const int nchains = b_start.n_cols;
	const bool free_calibration = all(fixed_b != fixed_b); // NA != NA
	const int max_cat = max(last-first)+1;
	
	progress_prl pb(ndraws, progress_init);	
	dqrng::xoshiro256plus rng(SEED);
	
	// cumulatives for bookkeeping
	const int nit = nbi.n_elem;
	ivec cnbi(nit+1);
	cnbi(0) = 0;
	cnbi.tail(nit) = cumsum(nbi);
	
	const int nbk = nib.n_elem;
	ivec cnib(nbk+1);
	cnib(0) = 0;
	cnib.tail(nbk) = cumsum(nib);
	
	ivec cbmax(nbk+1);
	cbmax[0] = 0;	
	for(int k=0; k<nbk; k++)
		cbmax[k+1] = cbmax[k] + bmax[k]+1;
	

	const int max_bscore = max(bmax);

	mat out_b(b_start.n_rows,ndraws), out_bklambda(bkscoretab.n_elem, ndraws);
	
#pragma omp parallel num_threads(max_cores)
	{
		const int thread = omp_get_thread_num();
		// working variables
		vec y(max_cat, fill::zeros), z(nbk, fill::zeros);
		vec bklambda(bkscoretab.n_elem, fill::ones);	// to do: fill 0 lijkt beter, onmogelijke scores?
		
		std::vector<long double> g(max_bscore+1);
		vec pi_k(max_bscore+1, fill::zeros);
	
		vec b(b_start.n_rows);

		int bk,out_index;
		long double sm;
		
		dqrng::xoshiro256plus lrng(rng);      		
		lrng.long_jump(thread + 1);	
		
#pragma omp for
		for(int chain=0; chain<nchains; chain++)
		{
			out_index = chain*(ndraws/nchains) + std::min(ndraws % nchains, chain);
			b = b_start.col(chain);
			const int niter = (ndraws/nchains + (chain < ndraws % nchains)-1)*step + warmup + 1;
			
			for (int iter=0; iter<niter; iter++)
			{
				if(pb.interrupted()) break;
				for (int k=0; k<nbk; k++)
				{
					elsym_helper(b,a, bfirst.memptr() + cnib[k], blast.memptr() + cnib[k], -1, nib[k], bmax[k], g);
					
					sm = 0;
					for (int s=0; s<=bmax[k];s++)
					{
					  if (g[s]>0){
						pi_k[s] = rgamma(lrng, bkscoretab[cbmax[k]+s]+0.1, 1.0);
						sm += pi_k[s];
					  }
					}
					
					for (int s=0; s<=bmax[k];s++)
					{
					  if(g[s]>0)
					  {
						pi_k[s] = pi_k[s]/sm;
						bklambda[cbmax[k]+s] = (pi_k[s]*m[k])/g[s];
					  }
					}
					
					z[k] = rgamma(lrng, m[k], 1.0/m[k]);
					
				}
				// bayes_items
				for (int i=0; i<nit;i++)
				{
					y.zeros();
					for (int ib_indx=cnbi[i]; ib_indx<cnbi[i+1]; ib_indx++)
					{
						bk = ib[ib_indx];
						elsym_helper(b, a, bfirst.memptr() + cnib[bk], blast.memptr() + cnib[bk], bi[ib_indx], nib[bk], bmax[bk], g);
							
						for (int s=0; s<=bmax[bk]-a[last[i]]; s++)
							for (int j=first[i]+1, c=0; j<=last[i]; j++,c++)
								y[c] += z[bk] * g[s] * bklambda[cbmax[bk]+s+a[j]];			
					}

					for (int j=first[i]+1, c=0; j<=last[i]; j++,c++)
						b[j] = rgamma(lrng, sufI[j] + prior_eta, 1/(y[c] + prior_rho));
				}
				
				if (!free_calibration)
				{
					b.elem(find_finite(fixed_b)) = fixed_b.elem(find_finite(fixed_b));
				}
				b.replace(datum::nan, 1); 		
				if(iter >= warmup && (iter - warmup) % step == 0)	
				{
					pb.tick(thread == 0);
					out_bklambda.col(out_index) = bklambda;
					out_b.col(out_index++) = b;
					if(thread==0) pb.checkInterrupt(); 
				}
			}
		}
	}	
	
	return Rcpp::List::create(Rcpp::Named("b") = out_b.t(), Rcpp::Named("lambda") = out_bklambda.t());
}

