#include<math.h>
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "elsymean.h"
#include "pascal.h"

// compute Symmetric Basis Functions without item (stored in gi)
// Warning: It is necessary that elements of a are montone nondecreasing !
// This is used in CML and should be replaced with time
void ElSym0(double *b, int *a, int *first, int *last, int *item1,int* item2, int *nI, int *mS, double *g)
{
  int Msc=0, initial=0;
  
  for (int s=0;s<=mS[0];s++) {g[s]=0;} // make sure gi is 0 to start with
  initial=0;
  if ((item1[0]==0)&&(item2[0]==-1)) {initial=1;}
  if ((item1[0]==-1)&&(item2[0]==0)) {initial=1;}
  if ((item1[0]>1)&&(item2[0]==0)) {initial=1;}
  if ((item1[0]==0)&&(item2[0]>1)) {initial=1;}
  if ((item1[0]==1)&&(item2[0]==0)) {initial=2;}
  if ((item1[0]==0)&&(item2[0]==1)) {initial=2;}
  
  g[0]=1;
  for (int j=first[initial];j<=last[initial];j++){g[a[j]]=b[j];}
  Msc=a[last[initial]];
  for (int i=1;i<nI[0];i++)
  {
    if ((!(i==item1[0]))&&(!(i==item2[0]))&&(!(i==initial)))
    {
      for (int s=Msc;s>=0;s--)
      {
        for (int j=last[i];j>=first[i];j--)
        {
          g[s+a[j]]+=g[s]*b[j];
        }
      }
      Msc+=a[last[i]];
    }
  }
}

//This elsym function allow b to contain the parameters that would normally be
//implicit and equal to one
void ElSym(double *b, int *a, int *first, int *last, int *item1, int* item2, int *nI, int *mS, double *g) 
{
  int Msc, idx;
 
  double *gg;
  gg=(double*)calloc((2*(mS[0]+1)), sizeof(double)); 
 /* double *gg=NULL;
  void *_gg = realloc(gg, (2*(mS[0]+1) * sizeof(double)));
  gg=(double*)_gg;*/
  
  idx=0; // start from column one
  gg[0]=1;
  Msc=0;
  for (int i=0;i<nI[0];i++)
  {
    if ((!(i==item1[0]))&&(!(i==item2[0])))
    {
      for (int s=Msc;s>=0;s--) {gg[2*s+(1-idx)]=0;}
      for (int s=Msc;s>=0;s--)
      {
        for (int j=first[i];j<=last[i];j++)
        {
          if (b[j]>0) { gg[2*(s+a[j])+(1-idx)]+=gg[2*s+idx]*b[j]; }
        }
      }
      Msc+=a[last[i]];
      idx=1-idx; // swap columns
    }
  }
  for (int s=0;s<=Msc;s++) 
  {
     g[s]=gg[2*s+idx];
  } 
  free(gg);
}



void meanElSym(double *b, int *a, int *first, int *last, int *item1,int* item2, int *nI, int *mS, double *g, int *idx)
{
  int Msc,n;
  double cij, gg[10000]={0};
  
  idx[0]=0; // start from column one
  gg[0]=1;
  Msc=0;
  for (int i=0;i<nI[0];i++)
  {
    if ((!(i==item1[0]))&&(!(i==item2[0])))
    {
      n=Msc+a[last[i]];
      for (int s=Msc;s>=0;s--) {gg[2*s+(1-idx[0])]=0;}
      for (int s=0;s<=Msc;s++)
      {
        for (int j=first[i];j<=last[i];j++)
        {
          if (b[j]>0)
          {
            cij=exp(lbinom(Msc,s)-lbinom(n,s+a[j]));
            gg[2*(s+a[j])+(1-idx[0])]+=gg[2*s+idx[0]]*b[j]*cij;
          }
        }
      }
      Msc=n;
      idx[0]=1-idx[0]; // swap columns
    }
  }
  for (int s=0;s<=Msc;s++) {g[s]=gg[2*s+idx[0]];}
}

void ittot_c(double *b, int *a, double *c, int *it, int *first, int *last, int *nI, int *ms, double *pi)
{
  int ncat    = last[it[0]]-first[it[0]]; // minus de laagste 
  int nscores = ms[0]+1;
  int npar    = last[nI[0]-1]+1;
  int ni, k, i1=-1, i2=-1, idx, indx;
  double lchoose_ns;
  
  double *g=NULL;
  void *_g = realloc(g, (nscores * sizeof(double)));
  g=(double*)_g;
  
  double *gi=NULL;
  void *_gi = realloc(gi, (nscores * sizeof(double)));
  gi=(double*)_gi;
  
  double *eta=NULL;
  void *_eta = realloc(eta, (npar * sizeof(double)));
  eta=(double*)_eta;
  
  double *logb=NULL;
  void *_logb = realloc(logb, npar * sizeof(double));
  logb=(double*)_logb;
  
  double *logc=NULL;
  void *_logc = realloc(logc, npar * sizeof(double));
  logc=(double*)_logc;
  
  // prepare log parameters
  for (int h=0;h<npar;h++)
  {
    logb[h]=log(b[h]);
    logc[h]=log(c[h]);
  }
  
  ni=ms[0]-a[last[it[0]]]; 
  for (int s=0;s<=ms[0];s++)
  {
    lchoose_ns = lbinom(ms[0],s);
    for (int h=0;h<npar;h++) { eta[h]=exp(logb[h]+(a[h]*s)*logc[h]); }
    meanElSym(eta, a, first, last, &i1, &i2, nI, ms, g, &idx);
    meanElSym(eta, a, first, last,  it, &i2, nI, ms, gi, &idx);
    k=0; 
    for (int j=(first[it[0]]+1);j<=last[it[0]]; j++) 
    {
      indx=s-a[j];
      if ( (indx>-1)&&(indx<(nscores-a[last[it[0]]])) )
      {
        pi[ncat*s+k]=log(eta[j])+log(gi[indx])-log(g[s]);
        pi[ncat*s+k]+= lbinom(ni,indx)-lchoose_ns; //lchoose(ni,indx)-lchoose(ms[0],s);
        pi[ncat*s+k]=exp(pi[ncat*s+k]);
      }
      k++;
    }
  }
  
  free(logb);
  free(logc); 
  free(g);
  free(gi);
  free(eta);
}

void ittot_mat0(double *b, int *a, double *c, int *first, int *last, int *nI, int *ms, double *pi)
{
  int nscores = ms[0]+1;
  int npar    = last[nI[0]-1]+1;
  int ni, k, i1=-1, i2=-1, indx;
  
  double *g=NULL;
  void *_g = realloc(g, (nscores * sizeof(double)));
  g=(double*)_g;
  
  double *gi=NULL;
  void *_gi = realloc(gi, (nscores * sizeof(double)));
  gi=(double*)_gi;
  
  double *eta=NULL;
  void *_eta = realloc(eta, (npar * sizeof(double)));
  eta=(double*)_eta;
  
  double *logb=NULL;
  void *_logb = realloc(logb, npar * sizeof(double));
  logb=(double*)_logb;
  
  double *logc=NULL;
  void *_logc = realloc(logc, npar * sizeof(double));
  logc=(double*)_logc;
  
  // prepare log parameters
  for (int h=0;h<npar;h++)
  {
    logb[h]=log(b[h]);
    logc[h]=log(c[h]);
  }
  
  for (int s=0;s<=ms[0];s++)
  {
    k=0; 
    for (int h=0;h<npar;h++) { eta[h]=exp(logb[h]+(a[h]*s)*logc[h]); }
    ElSym(eta, a, first, last, &i1, &i2, nI, ms, g);
    for (int it=0;it<nI[0];it++)
    {
      ElSym(eta, a, first, last,  &it, &i2, nI, ms, gi);
      ni=ms[0]-a[last[it]]; 
      for (int j=first[it];j<=last[it]; j++) 
      {
        indx=s-a[j];
        if ( (indx>=0)&&(indx<(nscores-a[last[it]])) )
        {
          pi[npar*s+k] =log(eta[j])+log(gi[indx])-log(g[s]);
          pi[npar*s+k] =exp(pi[npar*s+k]);
        }
        k++;
      }
    }
  }
  
  free(logb);
  free(logc); 
  free(g);
  free(gi);
  free(eta);
}


void ittot_mat(double *b, int *a, double *c, int *first, int *last, int *nI, int *ms, double *pi)
{
  int nscores = ms[0]+1;
  int npar    = last[nI[0]-1]+1;
  int ni, k, i1=-1, i2=-1, idx, indx;
  double lchoose_ns;
  
  double *g=NULL;
  void *_g = realloc(g, (nscores * sizeof(double)));
  g=(double*)_g;
  
  double *gi=NULL;
  void *_gi = realloc(gi, (nscores * sizeof(double)));
  gi=(double*)_gi;
  
  double *eta=NULL;
  void *_eta = realloc(eta, (npar * sizeof(double)));
  eta=(double*)_eta;
  
  double *logb=NULL;
  void *_logb = realloc(logb, npar * sizeof(double));
  logb=(double*)_logb;
  
  double *logc=NULL;
  void *_logc = realloc(logc, npar * sizeof(double));
  logc=(double*)_logc;
  
  // prepare log parameters
  for (int h=0;h<npar;h++)
  {
    logb[h]=log(b[h]);
    logc[h]=log(c[h]);
  }
  
  for (int s=0;s<=ms[0];s++)
  {
    k=0; 
    for (int h=0;h<npar;h++) { eta[h]=exp(logb[h]+(a[h]*s)*logc[h]); }
    lchoose_ns = lbinom(ms[0],s);
    meanElSym(eta, a, first, last, &i1, &i2, nI, ms, g, &idx);
    for (int it=0;it<nI[0];it++)
    {
      meanElSym(eta, a, first, last,  &it, &i2, nI, ms, gi, &idx);
      ni=ms[0]-a[last[it]]; 
      for (int j=first[it];j<=last[it]; j++) 
      {
        indx=s-a[j];
        if ( (indx>=0)&&(indx<(nscores-a[last[it]])) )
        {
          pi[npar*s+k] =log(eta[j])+log(gi[indx])-log(g[s]);
          pi[npar*s+k]+= lbinom(ni,indx)-lchoose_ns; 
          pi[npar*s+k] =exp(pi[npar*s+k]);
        }
        k++;
      }
    }
  }
  
  free(logb);
  free(logc); 
  free(g);
  free(gi);
  free(eta);
}





/// routines for CML //

void E_n(double *b, int *a, int *first, int *last, int *scoretab, int *nI, int *mS, double *expect)
{
  int i1,i2;
 /* int nscores = mS[0]+1;
  double *g=NULL;
  void *_g = realloc(g, (nscores * sizeof(double)));
  g=(double*)_g;
  double *gi=NULL;
  void *_gi = realloc(gi, (nscores * sizeof(double)));
  gi=(double*)_gi; */
 
  
  double g[10000]={0};
  double gi[10000]={0}; 
  
  // pre-compute Symmetric Basis Functions
  i1=-1;
  i2=-1;
  ElSym(b,a,first,last,&i1,&i2,nI,mS,g);
  for (int item=0;item<nI[0];item++)
  {
    ElSym(b,a,first,last,&item,&i2,nI,mS,gi);
    for (int j=first[item];j<=last[item];j++)
    {
      for (int s=a[j];s<=mS[0];s++)
      {
        if (g[s]>0) expect[j]+=scoretab[s]*gi[s-a[j]]*b[j]/g[s];
      }
    }
  }
/*  free(g);
  free(gi); */
}

void E(double *b, int *a, int *first, int *last, int *scoretab, int *nI, int *mS, double *expect)
{
  int i1,i2;
  double *g=NULL;
  double *gi=NULL;
  // initialize g and gi
  void *_g = realloc(g, ((mS[0]+1) * sizeof(double)));
  void *_gi = realloc(gi, ((mS[0]+1) * sizeof(double)));
  g = (double*)_g;
  gi = (double*)_gi;
  // pre-compute Symmetric Basis Functions
  i1=-1;
  i2=-1;
  ElSym0(b,a,first,last,&i1,&i2,nI,mS,g);
  for (int item=0;item<nI[0];item++)
  {
    ElSym0(b,a,first,last,&item,&i2,nI,mS,gi);
    for (int j=first[item];j<=last[item];j++)
    {
      for (int s=a[j];s<=mS[0];s++)
      {
        if (g[s]>0) expect[j]+=scoretab[s]*gi[s-a[j]]*b[j]/g[s];
      }
    }
  }
  free(g);
  free(gi);
}

/*
H: Routine to calculate the information matrix/Hessian
*/

void H(double *b, int *a, int* nPar, int *first, int *last, int *scoretab, int *nI, int *mS, double *Hess)
{
  int i1,i2,tel=-1, jump=0;
  double *g=NULL;
  double *gi=NULL;
  double *gk=NULL;
  double *gik=NULL;
  // initialize g and gi
  void *_g = realloc(g, ((mS[0]+1) * sizeof(double)));
  void *_gi = realloc(gi, ((mS[0]+1) * sizeof(double)));
  void *_gk = realloc(gk, ((mS[0]+1) * sizeof(double)));
  void *_gik = realloc(gik, ((mS[0]+1) * sizeof(double)));
  g = (double*)_g;
  gi = (double*)_gi;
  gk = (double*)_gk;
  gik = (double*)_gik;
  
  // pre-compute Symmetric Basis Functions
  i1=-1;i2=-1;
  ElSym0(b,a,first,last,&i1,&i2,nI,mS,g);
  // compute expected sufficient statistics
  for (int item=0;item<nI[0];item++)
  {
    tel+=1+jump;
    jump+=(last[item]-first[item]+1);
    
    // compute Symmetric Basis Functions without item (stored in gi)
    ElSym0(b,a,first,last,&item,&i2,nI,mS,gi);
    
    // compute expected sufficient statistics
    for (int j=first[item];j<=last[item];j++)
    {
      for (int s=a[j];s<=mS[0];s++)
      {
        if (g[s]>0)
        {
          Hess[j+nPar[0]*j]+=scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(1-(gi[s-a[j]]*b[j]/g[s]));
        }
      }
      tel++;
      for (int k=(j+1);k<=last[item];k++)
      {
        for (int s=a[k];s<=mS[0];s++)
        {
          if (g[s]>0)
          {
            Hess[k+nPar[0]*j]-=scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(gi[s-a[k]]*b[k]/g[s]);
          }
        }
        tel++;
      }
      for (int k=(item+1);k<nI[0];k++)
      {
        ElSym0(b,a,first,last,&k,&i2,nI,mS,gk);
        ElSym0(b,a,first,last,&k,&item,nI,mS,gik);
        for (int l=first[k];l<=last[k];l++)
        {
          for (int s=0;s<=mS[0];s++)
          {
            if (g[s]>0)
            {
              if (s>=(a[j]+a[l]))      { Hess[l+nPar[0]*j]+= scoretab[s]*(gik[s-a[j]-a[l]])*((b[j]*b[l])/g[s]);}
              if ((s>=a[j])&&(s>=a[l])) { Hess[l+nPar[0]*j]-=scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(gk[s-a[l]]*b[l]/g[s]);}
            }
          }
          tel++;
        }
      }
    }
  }
  free(g);
  free(gi);
  free(gk);
  free(gik);
}
