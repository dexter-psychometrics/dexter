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