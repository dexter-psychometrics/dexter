#include <RcppArmadillo.h>
#include "pascal.h"

// (versions of) routines that use meanElsym

// This calculates elsym functions divided by appropriate binomial coeficients using
// Pascal.c which stores the binomial coeficients
// Purpose is to keep the elsym functions from over- or underflowing
// This function also allows b to contain the parameters that would normally be
// implicit and equal to one


using namespace arma;


//[[Rcpp::export]]
void meanElSym(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, arma::vec& g) 
{
  int Msc=0, n, idx = 0; // start from column one
  const int nscores = mS + 1;
  double cij;
  vec gg(2 * nscores, fill::zeros);
  g.zeros();
  gg[0]=1;

  for (int i=0;i<nI;i++)
  {
    if (i != item1 && i != item2)
    {
      n=Msc+a[last[i]];
      for (int s=Msc;s>=0;s--) 
		gg[2*s+(1-idx)]=0;
      for (int s=0;s<=Msc;s++)
      {
        for (int j=first[i];j<=last[i];j++)
        {
          if (b[j]>0)
          {
            cij = exp(lbinom(Msc,s) - lbinom(n,s+a[j]));
            gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]*cij;
          }
        }
      }
      Msc=n;
      idx=1-idx; // swap columns
    }
  }
  for (int s=0;s<=Msc;s++) 
	g[s]=gg[2*s+idx];
}

void meanElSym12(const arma::vec& b, const arma::ivec& a, const arma::ivec& first, const arma::ivec& last, const int item1, const int item2, const int nI, const int mS, arma::vec& g1, arma::vec& g2) 
{
	int Msc=0, n, idx = 0; // start from column one
	const int nscores = mS + 1;
	double cij;
	vec gg(2 * nscores, fill::zeros);
	g1.zeros();
	g2.zeros();
	gg[0]=1;

	for (int i=0;i<nI;i++)
	{
		if (i != item1 && i != item2)
		{
			n=Msc+a[last[i]];
			for (int s=Msc;s>=0;s--) 
				gg[2*s+(1-idx)]=0;
			for (int s=0;s<=Msc;s++)
			 {
				for (int j=first[i];j<=last[i];j++)
				{
				  if (b[j]>0)
				  {
					cij = exp(lbinom(Msc,s) - lbinom(n,s+a[j]));
					gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]*cij;
				  }
				}
			}
			Msc=n;
			idx=1-idx; // swap columns
		}
	}
	// store g2  
	for (int s=0; s<=Msc; s++)
		g2[s] = gg[2*s + idx];
	
	// 1 loop remaining to add item2 to g1
	for (int s=Msc; s>=0; s--)
		gg[2*s+(1-idx)] = 0;
	
	n = Msc+a[last[item2]];
	
	for (int s=Msc; s>=0; s--)
		for (int j=first[item2]; j<=last[item2]; j++)
			if (b[j]>0)
			{
				cij = exp(lbinom(Msc,s) - lbinom(n,s+a[j]));
				gg[2*(s+a[j])+(1-idx)] += gg[2*s+idx]*b[j]*cij;
			}

	Msc += a[last[item2]];
	idx = 1-idx;
	
	// store g1  
	for (int s=0; s<=Msc; s++)
		g1[s] = gg[2*s + idx];
}



// expect should bet set to zero by caller
void Expect_mean(const vec& b, const ivec& a, int* ptr_first, int* ptr_last, int* ptr_scoretab, int nI, int n_score, vec& expect)
{
  const int i1=-1, i2=-1;
  const int mS = n_score - 1;
  const ivec first(ptr_first, nI, false, true), last(ptr_last, nI, false, true);
  const ivec scoretab(ptr_scoretab, n_score, false, true);	 
  
  vec g(n_score), gi(n_score);
  
  meanElSym(b,a,first,last,i1,i2,nI,mS,g);
  
  for (int item=0; item<nI; item++)
  {    
	meanElSym(b,a,first,last,item,i2,nI,mS,gi);
    for (int j=first[item];j<=last[item];j++)
      for (int s=a[j];s<=mS;s++)
        if (g[s]>0) 
			expect[j] += scoretab[s] * b[j] *  (gi[s-a[j]]/g[s]) * exp(lbinom(mS-a[last[item]],s-a[j]) - lbinom(mS,s)) ;
  }
}

// [[Rcpp::export]]
arma::vec E_booklets_mean(const arma::vec& b, const arma::ivec& a, 
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
		Expect_mean(b, a, first.memptr() + cnit(bl), last.memptr() + cnit(bl), scoretab.memptr() + cn_score(bl), nit(bl), n_score(bl), Et);
		E += Et;
	}
	return E;
}



void E_Hess_mean(const vec& b, const ivec& a, const ivec& first, const ivec& last, int* ptr_scoretab, int nI, int n_score, vec& expect, mat& Hess)
{

  
	const int i1=-1, i2=-1;
	const int nPar = b.n_elem;
	const int mS = n_score - 1;
	double tmp;
	
	vec gi_g(mS+1, fill::zeros); //cache for gi[s-a[j]] divided by g[s]
	
	const ivec scoretab(ptr_scoretab, n_score, false, true);	  
	
	vec g(mS+1, fill::zeros), gi(mS+1, fill::zeros), gk(mS+1, fill::zeros), gik(mS+1, fill::zeros);	
	
	meanElSym(b,a,first,last,i1,i2,nI,mS,g);

	for (int item=0;item<nI;item++)
	{
		meanElSym(b,a,first,last,item,i2,nI,mS,gi);
		for (int j=(first[item]+1);j<=last[item];j++)
		{
		  for (int s=a[j];s<=mS;s++)
		  {		    
			if (g[s]>0)
			{
			  tmp = b[j] * ( gi[s-a[j]] / g[s] ) * exp(lbinom(mS-a[last[item]],s-a[j]) - lbinom(mS,s));
			  gi_g[s] = scoretab[s] * tmp;		
			  expect[j] += gi_g[s]; // this line replaces function Expect
			  Hess.at(j,j) += gi_g[s] * (1-tmp); // within current item
			  
			}
		  }
		  for (int k=(j+1);k<=last[item];k++)
		  {
			for (int s=a[k];s<=mS;s++)
			{
			  if (g[s]>0)
			  {
				Hess.at(k,j) -= gi_g[s] * b[k] * (gi[s-a[k]]/g[s])  * exp(lbinom(mS-a[last[item]],s-a[k]) - lbinom(mS,s)); //within current item 
			  }
			}
		  }

		  for (int k=(item+1);k<nI;k++) //deze loopt over alle volgende items
		  {
			meanElSym12(b,a,first,last,k,item,nI,mS,gk,gik);
			
			for (int l=(first[k]+1);l<=last[k];l++)
			{
			  for (int s=0;s<=mS;s++)
			  {
				if (g[s]>0)
				{
				  if (s >= (a[j]+a[l]))
				  { 
					Hess.at(l,j) += scoretab[s] * b[j] * b[l] *  (gik[s-a[j]-a[l]] / g[s])  * exp(lbinom(mS-a[last[k]]-a[last[item]],s-a[j]-a[l]) - lbinom(mS,s)); //j is current item, l is andere(opvolgende) items
				  }
				  if ((s>=a[j])&&(s>=a[l])) 
				  { 
					Hess.at(l,j) -= gi_g[s] * b[l] * (gk[s-a[l]] / g[s])  * exp(lbinom(mS-a[last[k]],s-a[l]) - lbinom(mS,s));
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




// [[Rcpp::export]]
void NR_booklets_mean(const arma::vec& b, const arma::ivec& a, 
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
		E_Hess_mean(b, a, first, last, scoretab.memptr(), nit(0), n_score(0), EsufI, H);
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
					
			E_Hess_mean(bkb, bka, bkfirst, bklast, scoretab.memptr() + cn_score(bl), nit(bl), n_score(bl), Et, Ht);
		
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


// to~do: unit test, ittotmat0 equal to ittotmat, also for NR, expect
//[[Rcpp::export]]
arma::mat ittotmat_C(const arma::vec& b, const arma::ivec& a, const arma::vec& c, const arma::ivec& first, const arma::ivec& last, const arma::ivec& ps)
{
	const int ms = accu(a(conv_to<uvec>::from(last)));
	const int nI = last.n_elem;
	const int npar = accu(last-first) + nI;
	const int nscores = ms+1;
	const int i1=-1, i2=-1;  
	const vec logb = log(b);
	const vec alogc = a % log(c);  

	mat pi(npar,nscores,fill::zeros);
 
 
#pragma omp parallel
	{ 
		vec g(nscores), gi(nscores);
		vec eta(npar);
		int k, indx, ni;
		double lbinom_ms_s, log_g_s;
#pragma omp for	
		for (int s=0;s<=ms;s++)
		{
			if(ps[s] == 1)
			{
				k=0; 
				eta = exp(logb + s * alogc);
				lbinom_ms_s = lbinom(ms, s);		
				meanElSym(eta, a, first, last, i1, i2, nI, ms, g);
				log_g_s = log(g[s]);
				for (int it=0;it<nI;it++)
				{
				  meanElSym(eta, a, first, last,  it, i2, nI, ms, gi);
				  ni=ms-a[last[it]]; 
				  for (int j=first[it];j<=last[it]; j++) 
				  {
					indx=s-a[j];
					if ( (indx>=0)&&(indx<(nscores-a[last[it]])) )
					{
					  pi.at(k,s) = exp((log(eta[j]) + log(gi[indx]) - log_g_s) + (lbinom(ni,indx) - lbinom_ms_s));
					} 
					k++;
				  }
				}
			}
		}
	}
	return pi;
}



// score by score probability tables of two tests

// functions assumes A and B make up the complete set
// [[Rcpp::export]]
arma::mat ss_table_enorm_C(const arma::ivec& a, const arma::vec& b, 
						const arma::ivec& first, const arma::ivec& last,
						const arma::ivec& firstA, const arma::ivec& lastA, 
						const arma::ivec& firstB, const arma::ivec& lastB)
{
	const int mscA = accu(a.elem(conv_to<uvec>::from(lastA)));
	const int mscB = accu(a.elem(conv_to<uvec>::from(lastB)));
	const int Msc = mscA + mscB;	
	const int niA = firstA.n_elem, niB = firstB.n_elem;
	const int nI = niA + niB;
	int s_b;
	double lbinom_Msc_s;
	
	vec gA(mscA+1), gB(mscB+1), g(Msc+1);
	
	mat out(mscA+1, mscB+1);
	
	meanElSym(b,a,firstA,lastA,-1,-1, niA, mscA,gA);
    meanElSym(b,a,firstB,lastB,-1,-1, niB,mscB, gB);
    meanElSym(b,a,first,last,-1,-1, nI, Msc,  g);

    
	out.fill( -datum::inf);
	
    for (int s=0; s<=Msc;s++)
    {
      if (g[s]>0)
      {
		lbinom_Msc_s = lbinom(Msc, s);
        for (int s_a = std::max(0,s-mscB); s_a <= std::min(s,mscA); s_a++)
        {
          s_b = s - s_a;
          if (gA[s_a] > 0 && gB[s_b]>0)
          {
            out.at(s_a,s_b) = std::log(gA[s_a]) + std::log(gB[s_b]) - std::log(g[s]);
            out.at(s_a,s_b) += lbinom(mscA, s_a) + lbinom(mscB, s_b) - lbinom_Msc_s;
          }
        }
      }
    }
	return exp(out);
}



// [[Rcpp::export]]
arma::mat ss_table_im_C(const arma::ivec& a, const arma::vec& b, const arma::vec& c,
						const arma::ivec& first, const arma::ivec& last,
						const arma::ivec& firstA, const arma::ivec& lastA, 
						const arma::ivec& firstB, const arma::ivec& lastB)
{
	const int mscA = accu(a.elem(conv_to<uvec>::from(lastA)));
	const int mscB = accu(a.elem(conv_to<uvec>::from(lastB)));
	const int Msc = mscA + mscB;
	const int niA = firstA.n_elem, niB = firstB.n_elem;
	const int nI = niA + niB;
	
	int s_b;
	double lbinom_Msc_s;
	
	const vec logb = log(b), alogc = a % log(c);
	vec eta(b.n_elem);
	vec gA(mscA+1), gB(mscB+1), g(Msc+1);
	
	mat out(mscA+1, mscB+1);
	
	out.fill( -datum::inf);
		
    for (int s=0; s<=Msc; s++)
    {
		eta = exp(logb + s * alogc);		
		meanElSym(eta, a, first, last,-1,-1, nI, Msc,  g);
	  
		if (g[s]>0)		
		{
			meanElSym(eta, a, firstA, lastA,-1,-1, niA, mscA,  gA);
			meanElSym(eta, a, firstB, lastB,-1,-1, niB, mscB,  gB);		
			lbinom_Msc_s = lbinom(Msc, s);
			
			for (int s_a = std::max(0,s-mscB); s_a <= std::min(s,mscA); s_a++)
			{
				s_b = s-s_a;
				if (gA[s_a] > 0 && gB[s_b]>0)
				{
					out.at(s_a,s_b) = std::log(gA[s_a]) + std::log(gB[s_b]) - std::log(g[s]);
					out.at(s_a,s_b) += lbinom(mscA, s_a) + lbinom(mscB, s_b) - lbinom_Msc_s;
				}
			}
		}
    }

	return exp(out);
}



