#include<math.h>
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "enorm.h"

//This elsym function allow b to contain the parameters that would normally be
//implicit and equal to one

void ElSym(double *b, int *a, int *first, int *last, int *item1,int* item2, int *nI, int *mS, double *g, int *idx) {
  int Msc=0,Item_in;
  double gg[10000]={0};

  idx[0]=0; // start from column one
  gg[0]=1;
  Msc=0;
  for (int i=0;i<nI[0];i++)
  {
    Item_in=0;
    for (int j=first[i];j<=last[i];j++) if (b[j]!=0) Item_in++;
    if ((!(i==item1[0]))&&(!(i==item2[0]))&&(Item_in!=0))
    {
      for (int s=Msc;s>=0;s--) {gg[2*s+(1-idx[0])]=0;}
      for (int s=Msc;s>=0;s--)
      {
        for (int j=first[i];j<=last[i];j++)
        {
          if (b[j]>0) { gg[2*(s+a[j])+(1-idx[0])]+=gg[2*s+idx[0]]*b[j]; }
        }
      }
      Msc+=a[last[i]];
      idx[0]=1-idx[0]; // swap columns
    }
  }
  for (int s=0;s<=Msc;s++) {g[s]=gg[2*s+idx[0]];} }


