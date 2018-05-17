#include<math.h>
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "roger.h"


// sample the 0-based indices of responses to one item with (0-based) index i: 
// NOT adapted for use with b_0 and a_0 in
void sampleNRM0(double *theta, double *b, int *a, int *i, int *first, int *last, int *m, int *response)
{
  double *p=NULL, u;
  int k;
  GetRNGstate(); // get seed
  void *_p = realloc(p, ((last[i[0]]-first[i[0]]+2) * sizeof(double)));
  p=(double*)_p;
  for (int pers=0;pers<m[0];pers++)
  {
    p[0]=1;
    k=1;
    for (int j=first[i[0]];j<=last[i[0]];j++)
    {
      p[k]=p[k-1]+b[j]*exp(a[j]*theta[pers]);
      k++;
    }
    u=p[k-1]*runif(0,1);
    k=0;
    while (u>p[k]) {k++;}
    response[pers]=k;
  }
  PutRNGstate(); // put seed
  free(p);
}

// sample the 0-based indices of responses to one item with (0-based) index i: 
// adapted for use with b_0 and a_0 in
void sampleNRM(double *theta, double *b, int *a, int *i, int *first, int *last, int *m, int *response)
{
  double *p=NULL, u;
  int k;
  GetRNGstate(); // get seed
  void *_p = realloc(p, ((last[i[0]]-first[i[0]]+2) * sizeof(double)));
  p=(double*)_p;
  for (int pers=0;pers<m[0];pers++)
  {
    p[0]=b[first[i[0]]]*exp(a[first[i[0]]]*theta[pers]); 
    k=1;
    for (int j=first[i[0]]+1;j<=last[i[0]];j++)
    {
      p[k]=p[k-1]+b[j]*exp(a[j]*theta[pers]);
      k++;
    }
    u=p[k-1]*runif(0,1);
    k=0;
    while (u>p[k]) {k++;}
    response[pers]=k;
  }
  PutRNGstate(); // put seed
  free(p);
}

// sample test-scores: This routine is adapted for use with b_0 and a_0.
void sampleNRM2(double *theta, double *b, int *a, int *first, int *last, int *nI, int *m, int *score)
{
  double *p=NULL, u;
  int k=0,maxS=0;
  GetRNGstate(); // get seed
  for (int i=0;i<nI[0];i++) {if ((last[i]-first[i]+2)>maxS) {maxS=(last[i]-first[i]+2);}}
  void *_p = realloc(p, ((maxS+1) * sizeof(double)));
  p=(double*)_p;
  
  for (int pers=0;pers<m[0];pers++)
  {
    score[pers]=0;
    for (int i=0;i<nI[0];i++)
    {
      p[0]=b[first[i]]*exp(a[first[i]]*theta[pers]); 
      k=1;
      for (int j=first[i]+1;j<=last[i];j++) // note the +1
      {
        p[k]=p[k-1]+b[j]*exp(a[j]*theta[pers]);
        k++;
      }
      u=p[k-1]*runif(0,1);
      k=0;
      while (u>p[k]) {k++;}
      if (k>0) {score[pers]+=a[first[i]+k];}
    }
  }
  PutRNGstate(); // put seed
  free(p);
}

void PV0(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, double *theta)
{
  double atheta=0.0;
	int x;
	double *ptheta;
 	ptheta=&atheta;
	double *p=NULL, u;
	int maxS=0, k=0;
 	GetRNGstate(); // get seed
	for (int i=0;i<nI[0];i++) {if ((last[i]-first[i]+2)>maxS) {maxS=(last[i]-first[i]+2);}}
	void *_p = realloc(p, ((maxS+1) * sizeof(double)));
	p=(double*)_p;
	for (int person=0;person<nP[0];person++)
	{
		x=-1;
    while (x!=score[person])
		{
			atheta=rnorm(mu[pop[person]],sigma[pop[person]]);
			x=0;
			for (int i=0;i<nI[0];i++)
			{
				p[0]=1;
				k=1;
				for (int j=first[i];j<=last[i];j++)
				{
					p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
					k++;
				}
				u=p[k-1]*runif(0,1);
				k=0;
				while (u>p[k]) {k++;}
				if (k>0) {x+=a[first[i]+k-1];}
			}
		}
		theta[person]=atheta;
	}
  free(p);
 	PutRNGstate(); // put seed
}


// Adapted by timo for inclusion of zero category
/*
void PV1(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, double *theta)
{
  double atheta=0.0;
  int x;
  double *ptheta;
  ptheta=&atheta;
  double *p=NULL, u;
  int maxS=0, k=0;
  GetRNGstate(); // get seed
  for (int i=0;i<nI[0];i++) {if ((last[i]-first[i]+2)>maxS) {maxS=(last[i]-first[i]+2);}}
  void *_p = realloc(p, ((maxS+1) * sizeof(double)));
  p=(double*)_p;
  for (int person=0;person<nP[0];person++)
  {
    x=-1;
    while (x!=score[person])
    {
      atheta=rnorm(mu[pop[person]],sigma[pop[person]]);
      x=0;
      for (int i=0;i<nI[0];i++)
      {
        p[0]=b[first[i]]*exp(a[first[i]]*atheta); 
        k=1;
        for (int j=first[i]+1;j<=last[i];j++) // note the +1
        {
          p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
          k++;
        }
        u=p[k-1]*runif(0,1);
        k=0;
        while (u>p[k]) {k++;}
        if (k>0) {x+=a[first[i]+k];} // note that -1 has been deleted
      }
    }
    theta[person]=atheta;
  }
  free(p);
  PutRNGstate(); // put seed
}
 */

// Adapted by timo for inclusion of zero category
// produces nPV plausible values per score
void PV(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP, int *nI, int *nPop, int *nPV, double *theta)
{
  double atheta=0.0;
  int x;
  double *ptheta;
  ptheta=&atheta;
  double *p=NULL, u;
  int maxS=0, k=0;
  GetRNGstate(); // get seed
  for (int i=0;i<nI[0];i++) {if ((last[i]-first[i]+2)>maxS) {maxS=(last[i]-first[i]+2);}}
  void *_p = realloc(p, ((maxS+1) * sizeof(double)));
  p=(double*)_p;
  for (int t=0;t<nPV[0];t++)
  {
    for (int person=0;person<nP[0];person++)
    {
      x=-1;
      while (x!=score[person])
      {
        atheta=rnorm(mu[pop[person]],sigma[pop[person]]);
        x=0;
        for (int i=0;i<nI[0];i++)
        {
          p[0]=b[first[i]]*exp(a[first[i]]*atheta); 
          k=1;
          for (int j=first[i]+1;j<=last[i];j++) // note the +1
          {
            p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
            k++;
          }
          u=p[k-1]*runif(0,1);
          k=0;
          while (u>p[k]) {k++;}
          if (k>0) {x+=a[first[i]+k];} // note that -1 has been deleted
        }
      }
      theta[t*nP[0]+person]=atheta;
    }
  }
  free(p);
  PutRNGstate(); // put seed
}
/*
 Inputs n; a vector of length Mxs+1, each entry corresponding to a score with the initial value npv.
 If a score is not possible given the weights in a n = 0
 */
void recyclePV(double *b, int *a, int *first, int *last, int *nI, int *n, double *mu, double *sigma, double *R)
{
  int nzero=0, k=0, npv=n[0]; // zero can always occur
  double atheta, *p=NULL, u;
  int sm, ms, Mxs, colindx, nrows, iter;
  ms=0;Mxs=0;
  for (int i=0;i<nI[0];i++) {
    if ((last[i]-first[i]+2)>ms) {ms=(last[i]-first[i]+2);}
    Mxs += a[last[i]];
  }
  for (int i=1;i<=Mxs;i++){ 
    if (n[i]==0) nzero++;
  }
  
  void *_p = realloc(p, ((ms+1) * sizeof(double)));
  p=(double*)_p;
  GetRNGstate(); // get seed
  
  iter=0; nrows=Mxs+1;
  while (iter<((nrows-nzero)*npv))
  {
    atheta=rnorm(mu[0],sigma[0]);
    sm=0;
    for (int i=0;i<nI[0];i++)
    {
      p[0]=b[first[i]]*exp(a[first[i]]*atheta); 
      k=1;
      for (int j=first[i]+1;j<=last[i];j++) // note the +1
      {
        p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
        k++;
      }
      u=p[k-1]*runif(0,1);
      k=0;
      while (u>p[k]) {k++;}
      if (k>0) {sm+=a[first[i]+k];}
    }
   
    if (n[sm]>0) {
      colindx=npv-n[sm];
      R[sm+colindx*nrows]=atheta;
      n[sm]--;
      iter++;
    }
  }
  free(p);
  PutRNGstate(); // put seed
}

/*
 Inputs n; a vector of length Mxs+1, each entry corresponding to a score with value npv.
 This version accepts a vector called prior which contains a sample from the prior
*/
void recyclePV2(double *b, int *a, int *first, int *last, int *nI, int *n, double *prior, int *nprior, double *R)
{
  int nzero=0, k=0, npv=n[0], rindx;
  double atheta, *p=NULL, u;
  int sm, ms, Mxs, colindx, nrows, iter;
  ms=0;Mxs=0;
  for (int i=0;i<nI[0];i++) {
    if ((last[i]-first[i]+2)>ms) {ms=(last[i]-first[i]+2);}
    Mxs += a[last[i]];
  }
  for (int i=1;i<=Mxs;i++){ 
    if (n[i]==0) nzero++;
  }
  
  void *_p = realloc(p, ((ms+1) * sizeof(double)));
  p=(double*)_p;
  GetRNGstate(); // get seed
  
  iter=0;  nrows=Mxs+1;
  while (iter<((nrows-nzero)*npv))
  {
    rindx = (int)floor(runif(0,(nprior[0]-1))) ;
    atheta = prior[rindx];
    sm=0;
    for (int i=0;i<nI[0];i++)
    {
      p[0]=b[first[i]]*exp(a[first[i]]*atheta); 
      k=1;
      for (int j=first[i]+1;j<=last[i];j++) // note the +1
      {
        p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
        k++;
      }
      u=p[k-1]*runif(0,1);
      k=0;
      while (u>p[k]) {k++;}
      if (k>0) {sm+=a[first[i]+k];}
    }
    
    if (n[sm]>0) {
      colindx=npv-n[sm];
      R[sm+colindx*nrows]=atheta;
      n[sm]--;
      iter++;
    }
  }
  PutRNGstate(); // put seed
  free(p);
}

/*
 Inputs n; a vector of length Mxs+1, each entry corresponding to an A-weighted score with the initial value npv.
 If a score is not possible given the weights in a n = 0
 */
void recyclePVaA(double *b, int *a, int *A, int *first, int *last, int *nI, int *n, double *mu, double *sigma, double *R)
{
  int nzero=0, k=0, npv=n[0]; // zero can always occur
  double atheta, *p=NULL, u;
  int sm, ms, Mxs, colindx, nrows, iter;
  ms=0;Mxs=0;
  for (int i=0;i<nI[0];i++) {
    if ((last[i]-first[i]+2)>ms) {ms=(last[i]-first[i]+2);}
    Mxs += A[last[i]];
  }
  for (int i=1;i<=Mxs;i++){ 
    if (n[i]==0) nzero++;
  }
  
  void *_p = realloc(p, ((ms+1) * sizeof(double)));
  p=(double*)_p;
  GetRNGstate(); // get seed
  
  iter=0; nrows=Mxs+1;
  while (iter<((nrows-nzero)*npv))
  {
    atheta=rnorm(mu[0],sigma[0]);
    sm=0;
    for (int i=0;i<nI[0];i++)
    {
      p[0]=b[first[i]]*exp(a[first[i]]*atheta); 
      k=1;
      for (int j=first[i]+1;j<=last[i];j++) // note the +1
      {
        p[k]=p[k-1]+b[j]*exp(a[j]*atheta);
        k++;
      }
      u=p[k-1]*runif(0,1);
      k=0;
      while (u>p[k]) {k++;}
      if (k>0) {sm+=A[first[i]+k];}
    }
    
    if (n[sm]>0) {
      colindx=npv-n[sm];
      R[sm+colindx*nrows]=atheta;
      n[sm]--;
      iter++;
    }
  }
  free(p);
  PutRNGstate(); // put seed
}

///////////////////////////////
void Escore(double *theta, double *score, double *b, int *a, int *first, int *last, int *n)
{
  double denom, num;
  score[0]=0.0;
  for (int i=0; i<n[0]; i++)
  {
    num=0.0;
    denom=1.0;
    for (int j=first[i]+1;j<=last[i];j++) // note +1
    {
      num  +=a[j]*b[j]*exp(a[j]*theta[0]);
      denom+=     b[j]*exp(a[j]*theta[0]);
    }
    score[0]+=num/denom;
  }
}

