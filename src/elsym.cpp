#include <RcppArmadillo.h>



/*
* Notes
*
* omp.h and omp pragma's need to be encased in ifndef for Apple and other less able systems
*
* arma::vector is returned in R as a matrix with 1 column, need: [,1,drop=T] to make vectors again
*
* arguments of exported functions need the arma:: prefix otherwise not recognized in RcppExports
*
* get and put rngstate probably not necessary anymore, to do: check
*
* packagename_init.C is no longer necessary, automatically taken care of, use //[[Rcpp::export]]
* to make a function available in R
*
*/


using namespace arma;



// some magic to easily switch between arma::vec and std::vec in order to use std::vector<long double> and arma::vec somewhat interchangably
// to do: remove templates, always call elsym with long double, have version callable from r do a typecast
void fill_zeros(arma::vec& v){ v.zeros();}
void fill_zeros(std::vector<long double>& v){ std::fill(v.begin(), v.end(), 0); }

// for use when type depends on input
template <typename T>
T VecDbl(const int sz)
{
	T v(sz, 0);
	return v;
}

template <>
vec VecDbl<vec>(const int sz)
{
    vec v(sz, fill::zeros);
	return v;
}



// argument g can be either a std::vector<long double> or an arma::vec 
template <typename V_TYPE>
void ElSym(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS,  V_TYPE& g) 
{
  int Msc = 0, idx = 0; // start from column one
  
  V_TYPE gg = VecDbl<V_TYPE>(2 * mS + 2);
  fill_zeros(g);
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

// R cannot call templated functions
//[[Rcpp::export]]
void ElSym_C(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, arma::vec& g)
{
	ElSym(b,a,first,last,item1, item2,nI,mS,g);
} 

// comput the gamma functions g1(-item1) and g2(-item1, -item2)
template <typename V_TYPE>
void ElSym12(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, V_TYPE& g1, V_TYPE& g2) 
{
	int Msc = 0, idx = 0; // start from column one
	
	  
	V_TYPE gg = VecDbl<V_TYPE>(2 * mS + 2);
	fill_zeros(g1);
	fill_zeros(g2);
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



// to do: this one might need long doubles as well
// to do: possible scores in itt elsymmean
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
		vec g(nscores), gi(nscores);
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
							pi.at(k,s) = exp(log(eta[j]) + log(gi[indx]) - log(g[s]));
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
	  
	std::vector<long double> g(mS+1, 0), gi(mS+1, 0), gk(mS+1, 0), gik(mS+1, 0);
	
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

void elsym_helper(const vec& b, const ivec& a, int* ptr_first, int* ptr_last, const int item1, const int nI, const int mS, vec& g)
{
	const ivec first(ptr_first, nI, false, true), last(ptr_last, nI, false, true);
	
	ElSym(b, a, first, last, item1, -1, nI, mS, g);
}



// for now assumes lambda out false, can be added later
// fixed_b must be length b (i.e. what was b=NULL must now be b={NA,NA...}

// [[Rcpp::export]]
arma::mat calibrate_Bayes_C(const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, 
				 const arma::ivec& ib, const arma::ivec& bi, const arma::ivec& nbi, 
				 const arma::ivec& nib,
				 arma::ivec& bfirst, arma::ivec& blast,
				 const arma::ivec& bnscore, const arma::ivec& m,
				 const arma::ivec& sufI, const arma::ivec& bkscoretab,
				 const arma::vec& b_in,  const arma::vec& fixed_b, const int nIter)
{

	const bool free_calibration = all(fixed_b != fixed_b); // NA != NA

	// cumulatives for bookkeeping
	const int nI = nbi.n_elem;
	ivec cnbi(nI+1);
	cnbi(0) = 0;
	cnbi.tail(nI) = cumsum(nbi);
	
	const int nB = nib.n_elem;
	ivec cnib(nB+1);
	cnib(0) = 0;
	cnib.tail(nB) = cumsum(nib);
	
	ivec cbnscore(nB+1);
	cbnscore(0) = 0;	
	cbnscore.tail(nB) = cumsum(bnscore);

	const int max_bnscore = max(bnscore);

	// working variables
	vec y(sufI.n_elem, fill::zeros), z(nB, fill::zeros);
	vec bklambda(bkscoretab.n_elem, fill::ones);	
	
	double f, tmp;	
	int bkl;
	vec g(max_bnscore);
	vec fpwr(max_bnscore);
	fpwr[0] = 1;
	
	// copy b
	vec b(b_in.memptr(), b_in.n_elem);
	
	//output
	mat bx(b.n_elem, nIter);
	

	for (int iter=0; iter<nIter; iter++)
	{
		for (int bl=0; bl<nB; bl++)
		{
			
			// data augmentation
			elsym_helper(b,a, bfirst.memptr() + cnib[bl], blast.memptr() + cnib[bl], -1, nib[bl], bnscore[bl]-1, g);
			z[bl] = R::rgamma(m[bl], 1.0/accu(g.subvec(0, bnscore[bl]-1) % bklambda.subvec(cbnscore[bl], cbnscore[bl+1]-1)));
			
			// update lambda
			for(int i=0; i<bnscore[bl];i++)
				bklambda[cbnscore[bl] + i] = (g[i]>0) ? R::rgamma(bkscoretab[cbnscore[bl] + i] + 0.1, 1.0/(g[i]*z[bl])) : 0;

			// scale lambda such that g*lambda~scoretab
			tmp = accu(g.subvec(0, bnscore[bl]-1) % bklambda.subvec(cbnscore[bl], cbnscore[bl+1]-1));
			for(int i=0; i<bnscore[bl];i++)
				if(g[i]>0)
					bklambda[cbnscore[bl] + i] = m[bl] * bklambda[cbnscore[bl] + i]/tmp; 

			z[bl] = R::rgamma(m[bl], 1.0/accu(g.subvec(0, bnscore[bl]-1) % bklambda.subvec(cbnscore[bl], cbnscore[bl+1]-1)));
		}
		// bayes_items
		for(int i=0; i<nI;i++)
		{
			y.subvec(first[i], last[i]).zeros();		
			for(int ib_indx=cnbi[i]; ib_indx < cnbi[i+1]; ib_indx++)
			{
				bkl = ib[ib_indx];
				
				g.zeros();
				// need the number of the item in the present booklet				
				elsym_helper(b,a, bfirst.memptr() + cnib[bkl], blast.memptr() + cnib[bkl], bi[ib_indx], nib[bkl], bnscore[bkl]-1, g);
				
				//first category
				y[first[i]] += z[bkl] * accu(g.subvec(0, bnscore[bkl] - a[last[i]]-1) % bklambda.subvec(cbnscore[bkl], cbnscore[bkl+1] - a[last[i]]-1));
				// middle categories, if any
				for(int j=first[i]+1; j<last[i]; j++)
				{			
					y[j] += z[bkl] * accu(g.subvec(0, bnscore[bkl] - a[last[i]]-1) % bklambda.subvec(cbnscore[bkl]  + a[j], cbnscore[bkl+1] - a[last[i]] + a[j]-1));
				}			
				// last category (we assume each item has a zero score and at least one other score)
				y[last[i]] += z[bkl] * accu(g.subvec(0, bnscore[bkl] - a[last[i]]-1) % bklambda.subvec(cbnscore[bkl]  + a[last[i]], cbnscore[bkl+1]-1));		
			}
			
			for(int j=first[i]; j<=last[i]; j++)
				b[j] = R::rgamma(sufI[j] + 1.1, 1.0/y[j]);
		}

		// identify within items
		for (int i=0; i<nI; i++)
		{
			for(int j=first[i]+1; j<=last[i]; j++)
				b[j] /= b[first[i]];			
			b[first[i]]=1;
		}

    	// to do: differences in this last part don't seem to influence the outcome at all so not sure if correct, Timo to check the code
		if (free_calibration)
		{
			f = b[1]; 
			if(a[1] != 1)
				for(int i=0; i<nI; i++)
					for(int j = first[i]+1; j<= last[i]; j++)
						b[j] /= pow(f, ((double)(a[j]))/a[1]);

			// Lambda
			fpwr[1] = f;
			for(int s=2; s < max_bnscore; s++)
				fpwr[s] = fpwr[s-1] * f;
			
			for (int bl=0; bl<nB; bl++)
			{	
				bklambda.subvec(cbnscore[bl], cbnscore[bl+1]-1) %= fpwr.head(bnscore[bl]);
				
				// to do ask Timo, I assume >0 works, R version has != 0 but that is a bit problematic in C
				if (bklambda[cbnscore[bl]] > 0)
					bklambda.subvec(cbnscore[bl], cbnscore[bl+1]-1) /= bklambda[cbnscore[bl]];
			}
		}
		else
		{
			b.elem(find_finite(fixed_b)) = fixed_b.elem(find_finite(fixed_b));
		}
		b.replace(datum::nan, 1); // deal with items that are not in data (this is currently impossible I think)
		bx.col(iter) = b;
	}
	return bx.t();
}



/* *********************

experimental
expects first and last to be 0-indexed

Timo's H_IM

enjoy!

*********************** */

// possible scores omitted for the moment
template <bool diagonal>
void H_im_tmpl(const arma::ivec& a, const arma::vec& b, const arma::vec& c, const arma::ivec& first, const arma::ivec& last, const arma::ivec& sufI, const arma::ivec& sufC, const arma::ivec& scoretab,
			arma::mat& H, arma::vec& Grad, arma::mat& pi_s)
{
	const int nI = first.n_elem;
	const ivec nCat = last-first+1;	
	const int nscores = scoretab.n_elem;
	const int ms = nscores-1;
  
	const vec logb = log(b);
	const vec alogc = a % log(c);
	
	for (int i=0; i<nI; i++)
	{
		Grad[last[i]] = sufC[i];
		for (int j = first[i]+1; j<=last[i]; j++)
		{
		  Grad[j-1] = sufI[j];
		}
	}

#pragma omp parallel
	{
		std::vector<long double> g(nscores), gi(nscores), gik(nscores);
		vec eta(b.n_elem);
#pragma omp for
		for (int s=0; s<=ms; s++)
		{
			eta = exp(logb + s*alogc);
			ElSym(eta,a,first,last,-1,-1,nI,ms,g);

			for (int i=0; i<nI; i++)
			{
				ElSym(eta,a,first,last,i,-1,nI,ms,gi);

				for (int j = first[i]; j<=last[i]; j++)
				{
					int idx = s - a[j];
					if (idx>=0 && idx <= ms-a[j])
					{
						pi_s.at(j,s) = exp(log(eta[j]) + log(gi[idx]) - log(g[s]));
					}
				}
			}
	  
			for (int i=0; i<nI; i++)
			{
				for (int k=0; k<nI; k++)
				{
					if (k==i)
					{
						double sm_j1 = 0.0, sm_j2 = 0.0;
						for (int j=first[i]+1; j<=last[i]; j++)
						{
							sm_j1 += a[j]*a[j]*pi_s.at(j,s);
							sm_j2 += a[j]*pi_s.at(j,s);
							double sm_h = 0;
							for (int h=first[k]+1; h<=last[k]; h++)
							{
								sm_h += a[h]*(pi_s.at(j,s)*pi_s.at(h,s));
								if (j==h)
								{
									H.at(j-1,h-1) +=  scoretab[s]*(pi_s.at(j,s)*(1-pi_s.at(j,s)));
								}
								else
								{
									H.at(j-1,h-1) -= scoretab[s]*(pi_s.at(j,s)*pi_s.at(h,s));
								}
							}
							H.at(j-1,last[i]) += s * scoretab[s]*(a[j]*pi_s.at(j,s) - sm_h);
							H.at(last[i],j-1) = H.at(j-1,last[i]);
						}
						H.at(last[i],last[i]) += s * s * scoretab[s] * (sm_j1 - sm_j2 * sm_j2);
					}
					else if(!diagonal)
					{
						ElSym(b,a,first,last,i,k,nI,ms,gik);
						double sm_jh = 0;
						int ms_ik = ms - a[last[i]] - a[last[k]];
						for (int j=first[i]+1; j<=last[i]; j++)
						{
							double sm_h = 0;
							for (int h=first[k]+1; h<=last[k]; h++)
							{
								double pi_ijhk = 0;
								int idx = s - a[j] - a[h];
								if (idx>=0 && idx <= ms_ik)
								{
									if (g[s] != 0) 
										pi_ijhk = (gik[idx]*eta[j]*eta[h]) / g[s];
									sm_jh += a[j] * a[h] * (pi_ijhk - pi_s.at(j,s) * pi_s.at(h,s));
									sm_h += a[h] * (pi_ijhk - pi_s.at(j,s) * pi_s.at(h,s));
								}
								// delta_i delta_k
								H.at(j-1,h-1) -= scoretab[s] * (pi_ijhk - pi_s.at(j,s)*pi_s.at(h,s));
							}
							H.at(last[k],j-1) -= s * scoretab[s] * sm_h; // \delta-\sigma
							H.at(j-1,last[k]) = H.at(last[k],j-1);
						}
						H.at(last[i],last[k]) -= s * s * scoretab[s] * sm_jh; // sigma/sigma
					}
				}
			}
		}
	}

	// compute gradient sequentially
	for (int s=0; s<=ms; s++)
		for (int i=0; i<nI; i++)
		{
			double sm_j2 = 0;
			for (int j = first[i]+1; j<=last[i]; j++)
			{
				Grad[j-1] = Grad[j-1] - scoretab[s] * pi_s.at(j,s); 
				sm_j2 = sm_j2 + a[j]*pi_s.at(j,s);
			}			
			Grad[last[i]] = Grad[last[i]] - s * scoretab[s] * sm_j2;
		}
	
	
    // identification
	if (!diagonal)
	{
		H.col(0).zeros();
		H.row(0).zeros();
		Grad[0] = 0;
		H.at(0,0) = 1;
	}
}


/*
	Cpp:
	mat H(b.n_elem, b.n_elem, fill::zeros);
	vec Grad(b.n_elem, fill::zeros);
	mat pi_s(b.n_elem, scoretab.n_elem);
	
	R:
	H = matrix(0, length(b), length(b))
	Grad = double(length(b))
	pi_s = matrix(0, length(b), length(scoretab))	

	for the moment without use of possible scores
	c is expected to have length equal to b
*/

// [[Rcpp::export]]
void H_im(const arma::ivec& a, const arma::vec& b, const arma::vec& c, const arma::ivec& first, const arma::ivec& last, const arma::ivec& sufI, const arma::ivec& sufC, const arma::ivec& scoretab,
			arma::mat& H, arma::vec& Grad, arma::mat& pi_s,
			const bool diagonal=false)
{	
	if(diagonal)
	{
		H_im_tmpl<true>(a, b, c, first, last, sufI, sufC, scoretab, H, Grad, pi_s);
	} 
	else
	{
		H_im_tmpl<false>(a, b, c, first, last, sufI, sufC, scoretab, H, Grad, pi_s);
	}
}






