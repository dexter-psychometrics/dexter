### Native R versions of C functions in Dexter ###

swap <- function(a)
{
  a[1] = a[1] + a[2]
  a[2] = a[1] - a[2]
  a[1] = a[1] - a[2]
  return(a)
}

elsymR<-function(b,a,first,last)
{
  n=length(first)
  Msc=0
  gg=matrix(0,sum(a[last])+1,2)
  gg[1,1]=1
  col=c(1,2) 
  for (i in 1:n)
  {
    for (s in 0:Msc) gg[s+1,col[2]]=0
    for (s in 0:Msc)
    {
      for (j in first[i]:last[i])
      {
        gg[s+1+a[j],col[2]]=gg[s+1+a[j],col[2]]+gg[s+1,col[1]]*b[j]
      }
    }
    Msc = Msc+a[last[i]];
    col=swap(col)
  }
  g=NULL
  for (s in 0:Msc) g=c(g,gg[s+1,col[1]]) 
  return(g)
}

## original in R
## Using elsym
ittotmat0R = function (b, c, a, first, last) 
{
  ms=sum(a[last])
  mm=sum(last-first+1)
  pi=matrix(0, mm,(ms+1))
  logb = log(b)
  logc = log(c)
  for (s in 0:ms)
  {
    eta = exp(logb + (a * s) * logc)
    g = elsym(eta, a, first, last)
    k = 1
    for (it in 1:length(first)) 
    {
      gi = elsym(eta, a, first[-it], last[-it])
      for (j in first[it]:last[it]) 
      {
        idx = s + 1 - a[j]
        if ((idx > 0) & (idx <= length(gi))) 
        {
          pi[k, s + 1] = exp(log(eta[j]) + log(gi[idx]) - log(g[s + 1]))
        }
        k = k + 1
      }
    }
  }
  return(pi)
}


recycle_pvR = function(b, a, first, last, scores=NULL, npv=1, mu=0, sigma=2)
{
  ms.a=sum(a[last])
  if (is.null(scores)) {
    scores=0:ms.a
  }else
  {
    if (!identical(intersect(0:ms.a,scores),scores)) stop("wrong scores in recycle_pvR")
  }
  n=rep(npv,length(scores))
  R=matrix(0,length(scores),npv)
  while (any(n>0))
  {
    atheta=rnorm(1,mu,sigma)
    sm=rscore(b,a,atheta,first,last)+1
    if (n[sm]>0)
    {
      R[sm,n[sm]]=atheta
      n[sm]=n[sm]-1
    }
  }
  return(R)
}

#expected test score given theta
# assumes 0 category included
EscoreR <- function(theta, b, a, first, last)
{
  n=length(first)
  escore=0
  for (i in 1:n)
  {
    num=0
    denom=1
    for (j in (first[i]+1):last[i])
    {
      num = num   + a[j]*b[j]*exp(a[j]*theta)
      denom = denom +      b[j]*exp(a[j]*theta)
    }
    escore = escore + num/denom
  }
  return(escore)
}

H.STEP_R<-function(b,a,first,last, scoretab)
{
  nI=length(last)
  ms=length(scoretab)-1
  H=matrix(0,length(a),length(a))
  
  g=elsym(b,a,first,last)
  for (item in 1:nI)
  {
    gi=elsym(b,a,first[-item],last[-item])
    for (j in (first[item]+1):last[item])
    {
      for (s in (a[j]+1):(ms))
      {
        if (g[s]>0)
        {
          H[j,j] = H[j,j]+scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(1-(gi[s-a[j]]*b[j]/g[s]))
        }
      }
      
      if ((j+1)<=last[item])
      {
        for (k in (j+1):last[item])
        {
          for (s in (a[k]+1):(ms))
          {
            if (g[s]>0)
            {
              H[k,j] = H[k,j]-scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(gi[s-a[k]]*b[k]/g[s]);
            }
          }
        }
      }
      
      if ((item+1)<=nI)
      {
        for (k in (item+1):nI)
        {
          gk=elsym(b,a,first[-k],last[-k])
          gik=elsym(b,a,first[-c(item,k)],last[-c(item,k)])
          for (l in (first[k]+1):last[k])
          {
            for (s in 1:ms)
            {
              if (g[s]>0)
              {
                if ((s>(a[j]+a[l]))&((s-a[j]-a[l])<=length(gik))){
                  H[l,j] = H[l,j] + scoretab[s]*(gik[s-a[j]-a[l]])*((b[j]*b[l])/g[s]) 
                }
                if ((s>a[j])&(s>a[l])) {
                  H[l,j] = H[l,j] - scoretab[s]*(gi[s-a[j]]*b[j]/g[s])*(gk[s-a[l]]*b[l]/g[s])
                }
              }
            }
          }
        }
        
      }
    }
  }
  H=H+t(H)
  diag(H)=diag(H)/2
  return(H)
}
