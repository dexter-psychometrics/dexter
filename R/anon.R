##################### Function relating to ability estimation
## Plausible Values
# pv1 <- function(b,a,first,last,score,mu=0,sigma=2)
# {
#   tmp=.C("PV1",
#          as.double(b),
#          as.integer(a),
#          as.integer(first-1),
#          as.integer(last-1),
#          as.double(mu),
#          as.double(sigma),
#          as.integer(score),
#          as.integer(rep(0,length(score))),
#          as.integer(length(score)),
#          as.integer(length(first)),
#          as.integer(1),
#          as.double(0*score))[[12]]
#   return(tmp)
# }

#void PV(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, int *nPV, double *theta);
pv <- function(b,a,first,last,score,npv,mu=0,sigma=2)
{
  tmp=.C("PV",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(mu),
         as.double(sigma),
         as.integer(score),
         as.integer(rep(0,length(score))),
         as.integer(length(score)),
         as.integer(length(first)),
         as.integer(1),
         as.integer(npv),
         as.double(rep(0*score,npv)))[[13]]
  tmp=as.vector(tmp)
  dim(tmp)=c(length(score),npv)
  return(tmp)
}

# simulate responses to a single item. Adapted for inclusion zero
renorm <- function(x,b,a,theta,first,last,i)
{
  m=length(theta)
  tmp=.C("sampleNRM",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(i-1),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(m),
         as.integer(x))[[8]]
  return(tmp)
}
# simulate responses to a single item. NOt adapted for inclusion of parameter for 
# zero category
renorm0 <- function(x,b,a,theta,first,last,i)
{
  m=length(theta)
  tmp=.C("sampleNRM0",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(i-1),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(m),
         as.integer(x))[[8]]
  return(tmp)
}

#void Escore(double *theta, double *score, double *b, int *a, int *first, int *last, int *n);
E_score=function(theta,b,a,first,last)
{
  escore=as.double(0.0)
  n=length(first)
  tmp=.C("Escore",
         as.double(theta),
         as.double(escore),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(n))[[2]]
  return(tmp)
}

# computes distribution of score weighted with A conditionally on score weighted with a
# used for estimation of theta from the unweighted score (for instance)
# colSums(G[,,(1-idx)+1]) are elementary symmetric function with a
elsymat <- function(b,a,A,first,last)
{
  ms.a=sum(a[last])
  ms.A=sum(A[last])
  n=length(first)
  G=array(0, c(ms.A+1, ms.a+1, 2))
  # init
  G[1,1,1]=1
  for (j in first[1]:last[1])
  {
    G[A[j]+1,a[j]+1,1]=b[j]
  }
  # recursion
  idx=1
  ms.a=a[last[1]]
  ms.A=A[last[1]]
  for (i in 2:n)
  {
    G[,,idx+1]=G[,,(1-idx)+1]
    for (c in 1:(ms.a+1))
    {
      for (r in 1:(ms.A+1))
      {
        #G[r,c,idx+1]=G[r,c,(1-idx)+1]
        for (j in first[i]:last[i])
        {
          G[r+A[j],c+a[j],idx+1]=G[r+A[j],c+a[j],idx+1]+G[r,c,(1-idx)+1]*b[j]
        }
      }
    }
    idx=1-idx
    ms.a=ms.a+a[last[i]]
    ms.A=ms.A+A[last[i]]
  }
  return(G[,,(1-idx)+1])
}

# ML estimation of theta
theta_MLE <- function(b,a,first,last)
{
  ms.a=sum(a[last])
  theta=rep(0,ms.a-1)
  for (s in 1:(ms.a-1))
  {
    escore=-1
    while (abs(escore-s)>1e-1)
    {
      escore=E_score(theta[s],b,a,first,last)
      theta[s]=theta[s]+log(s/escore)
    }
  }
  return(c(-Inf,theta,Inf))
}

# EAP based on npv plausible values
theta_EAP <- function(b,a,first,last,score,npv=20,mu=0,sigma=2)
{
  PV=pv(b,a,first,last,score,npv,mu=0,sigma=2)
  return(rowMeans(PV))
}

# ML estimation (using EM) of theta based on A
theta_aA <- function(b,a,A,first,last)
{
  n=length(first)
  ms.a=sum(a[last])
  ms.A=sum(A[last])
  theta=rep(0,ms.A-1)
  G=elsymat(b,a,A,first,last)
  g=colSums(G)
  s_set=which(rowSums(G)>0)-1
  s_set=s_set[-c(1,length(s_set))]
  for (s in s_set)
  {
    converged=FALSE
    while (!converged)
    {
      # E-step (expected score weighted with a, given score weighted with A and current theta)
      p=g*exp((0:ms.a)*theta[s])
      p=p/sum(p)
      E=(G[s+1,]/g)*p
      E=sum((0:ms.a)*E/sum(E))
      # M-step
      converged=(abs(log(E/E_score(theta[s],b,a,first,last)))<1e-6)
      theta[s]=theta[s]+log(E/E_score(theta[s],b,a,first,last))
    }
  }
  return(c(-Inf,theta,Inf))
}

# Expected distribution given one ability theta
pscore <- function(theta,b,a,first,last)
{
  g=elsym(b,a,first,last)
  p=rep(0.0,length(g))
  for (s in 1:length(g))
  {
    p[s]=g[s]*exp((s-1)*theta)
  }
  return(p/sum(p))
}

# estimate a single ability for a whole score distribution
# to be used for the 3DC standard setting method
# and testing for overdispersion
theta_score_distribution <- function(b,a,first,last,scoretab)
{
  ms.a=sum(a[last])
  theta=0
  np=sum(scoretab)
  escore=-1
  score=(0:ms.a)%*%scoretab
  while (abs(escore-score)>1e-6)
  {
    escore=np*E_score(theta,b,a,first,last)
    theta=theta+log(score/escore)
  }
  return(theta)
}


##################################
## Score-by-score table. Currently using mean_ElSym
## as with estim.. call with c=ic[C]
SSTable <- function(m, AB, model) {
  if (model=="IM") {ic=m$est$cIM; b=m$est$bIM} else {ic=m$est$cRM; b=m$est$bRM}
  first = m$ss$il$first
  last =  m$ss$il$last
  C = rep(1:nrow(m$ss$il), m$ss$il$nCat)
  a = m$ss$sl$item_score
  ic = ic[C]
  A = AB[[1]]
  B = AB[[2]]
  ### Check
  if (length(intersect(A,B))!=0) stop("sets not disjunct")

  ## Bookkeeping
  bA=NULL; bB=NULL
  aA=NULL; aB=NULL
  cA=NULL; cB=NULL
  lastA=NULL; firstA=NULL
  lastB=NULL; firstB=NULL
  telAF=1; telBF=1
  for (i in 1:length(first))
  {
    if (i %in% A)
    {
      bA=c(bA,b[first[i]:last[i]])
      aA=c(aA,a[first[i]:last[i]])
      cA=c(cA,ic[first[i]:last[i]])
      firstA=c(firstA,telAF)
      lastA=c(lastA,telAF+last[i]-first[i])
      telAF=telAF+last[i]-first[i]+1
    }
    if (i %in% B)
    {
      bB=c(bB,b[first[i]:last[i]])
      aB=c(aB,a[first[i]:last[i]])
      cB=c(cB,ic[first[i]:last[i]])
      firstB=c(firstB,telBF)
      lastB=c(lastB,telBF+last[i]-first[i])
      telBF=telBF+last[i]-first[i]+1
    }
  }
  MscA=sum(aA[lastA])
  MscB=sum(aB[lastB])
  Msc=sum(a[last])
  out=matrix(NA,MscA+1,MscB+1)
  
  ### Rasch Model
  if (model=="RM")
  {
    gA = mean_ElSym(bA,aA,firstA,lastA)
    gB = mean_ElSym(bB,aB,firstB,lastB)
    g =  mean_ElSym(b,a,first,last)
    
    for (s in 0:Msc)
    {
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b=s-s_a
        out[s_a+1,s_b+1] = log(gA[s_a+1]) + log(gB[s_b+1]) - log(g[s+1])
        out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
      }
    }
  }
  
  ### IM
  if (model=="IM")
  {
    logb=log(b); logc=log(ic)
    logbA=log(bA); logcA=log(cA)
    logbB=log(bB); logcB=log(cB)
    for (s in 0:Msc)
    {
      eta=exp(logb+(a*s)*logc)
      etaA=exp(logbA+(aA*s)*logcA)
      etaB=exp(logbB+(aB*s)*logcB)
      gA = mean_ElSym(etaA,aA,firstA,lastA)
      gB = mean_ElSym(etaB,aB,firstB,lastB)
      g =  mean_ElSym(eta,a,first,last)
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b=s-s_a
        out[s_a+1,s_b+1]=log(gA[s_a+1])+log(gB[s_b+1])-log(g[s+1])
        out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
      }
    }
  }
  return(list(tbl=exp(out),m=m,AB=AB,model=model))
}


#################################### Item-Total Regressions
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

## Wrapper to C function
# using Elsym
ittotmat0 <- function(b,c,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  mm=sum(last-first+1)
  pi=double(mm*(ms+1))
  tmp=.C("ittot_mat0",
         as.double(b),
         as.integer(a),
         as.double(c),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(n),
         as.integer(ms),
         as.double(pi))
  dim(tmp[[8]])=c(mm,(ms+1))
  return(as.matrix(tmp[[8]])) 
}

## Wrapper to C function
# currently using mean_ElSym
ittotmat <- function(b,c,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  mm=sum(last-first+1)
  pi=double(mm*(ms+1))
  tmp=.C("ittot_mat",
          as.double(b),
          as.integer(a),
          as.double(c),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(n),
          as.integer(ms),
          as.double(pi))
  dim(tmp[[8]])=c(mm,(ms+1))
  return(as.matrix(tmp[[8]])) 
}


###################################
# Currently with elsym and scale
EstIM  <- function(ss) {
  first = ss$il$first
  last = ss$il$last
  a = ss$sl$item_score
  sufI = ss$sl$sufI
  sufC = ss$il$sufC
  C = rep(1:nrow(ss$il), ss$il$nCat)
  scoretab = ss$tl$N
  
  m=sum(scoretab) ##
  nI=length(last)
  b=rep(0,length(sufI))
  ic=rep(1,length(sufC))
  se.ic=vector("numeric", nI)
  HRM=matrix(0,length(b),length(b))
  
  # Identification
  b[sufI>0]=1
  upd_set=vector("list",nI)
  for (i in 1:nI)
  {
    upd_set[[i]]=which(sufI[first[i]:last[i]]>0)
    upd_set[[i]]=upd_set[[i]][-1]
    upd_set[[i]]=(first[i]:last[i])[upd_set[[i]]]
  }
  
  converged=2
  while(converged>0.001)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last) ##
    for (i in 1:nI)
    {
      if (length(upd_set[[i]])>0)
      {
        pi=pi_mat[upd_set[[i]],,drop=FALSE]
        pi[is.na(pi)]=0
        E=sufI[upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # NR update for parameters of item i
        update=solve(H,E)
        b[upd_set[[i]]]=b[upd_set[[i]]]*exp(update)
        converged=max(converged,max(abs(E))/m) #
        HRM[upd_set[[i]],upd_set[[i]]]=H
      }
    }
  }
  
  bRM=b
  cRM=ic
  
  ## IM
  converged=2
  scale=2
  while(converged>0.001)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last) ##
    for (i in 1:nI)
    {
      # gradient and hessian for thresholds of item i
      if (length(upd_set[[i]])>0)
      {
        pi=pi_mat[upd_set[[i]],,drop=FALSE]
        pi[is.na(pi)]=0
        E=sufI[upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # gradient and hessian for interaction parameter
        ncol_pi=ncol(pi)
        E=c(E,sufC[i])
        H=cbind(H,rep(0,nrow(H)))
        H=rbind(H,rep(0,ncol(H)))
        k=1
        e0=0; e1=0
        f=matrix(0,nrow(pi),ncol_pi)
        g=matrix(0,nrow(pi),ncol_pi)
        h=0
        for (j in upd_set[[i]])
        {
          E[length(E)]=E[length(E)]-a[j]*sum((0:(ncol_pi-1))*scoretab*pi[k,])
          e0=e0+a[j]*pi[k,]
          e1=e1+a[j]^2*pi[k,]
          f[k,]=a[j]*(0:(ncol_pi-1))*pi[k,]
          g[k,]=pi[k,]
          h=h+a[j]*(0:(ncol_pi-1))*pi[k,]
          k=k+1
        }
        H[nrow(H),nrow(H)]=sum((0:(ncol_pi-1))^2*(e1-e0^2)*scoretab)
        for (k in 1:nrow(f))
        {
          H[k,nrow(H)]=sum((f[k,]-g[k,]*h)*scoretab)
          H[nrow(H),k]=H[k,nrow(H)]
        }
        # NR update for parameters of item i
        update=solve(H*scale,E)
        b[upd_set[[i]]]=b[upd_set[[i]]]*exp(update[-length(update)])
        ic[i]=ic[i]*exp(update[length(update)])
        se.ic[i]=solve(H)[nrow(H),nrow(H)]
        converged=max(converged,max(abs(E))/m)
      }
    }
    if (converged<1) scale=1
  }
  
  return(list(group=ss$group,bRM=bRM,cRM=cRM,bIM=b,cIM=ic,se.c=se.ic,HRM=HRM))
}


##################################################### calibrate

#### Bayes
calibrate = function (itemList, booklet, sufI, b, a, first, last, nIter) {
  nb = length(booklet)
  n = length(itemList)
  y = rep(0, length(sufI))
  z = NULL
  bx = matrix(0, nIter, length(b))
  lx=list()
  length(lx)=nb
  names(lx)=names(booklet)
  for (bl in 1:nb) lx[[bl]]=matrix(0,nIter,sum(a[booklet[[bl]]$last])+1)
  pb = txtProgressBar(min=0, max=nIter)
  for (iter in 1:nIter)
  {
    for (bl in 1:nb)
    {
      # data augmentation
      g = mean_ElSym(b, a, booklet[[bl]]$first, booklet[[bl]]$last)
      lg1 = length(g) - 1
      scale_g=choose(lg1, 0:lg1)
      z[bl] = rgamma(1, shape=booklet[[bl]]$m, rate=sum(g*scale_g*booklet[[bl]]$lambda))
      # update lambda
      idx = which(g != 0.0)
      booklet[[bl]]$lambda[idx] = rgamma(length(idx), shape=booklet[[bl]]$scoretab[idx]+1.1, rate=(g*scale_g*z[bl])[idx])
      booklet[[bl]]$lambda[-idx] = 0.0
      # scale lambda such that g*lambda~scoretab
      booklet[[bl]]$lambda[idx] = booklet[[bl]]$m*booklet[[bl]]$lambda[idx]/sum(g*scale_g*booklet[[bl]]$lambda)
      z[bl] = rgamma(1, shape=booklet[[bl]]$m, rate=sum(g*scale_g*booklet[[bl]]$lambda))
    }
    for (i in 1:n)
    {
      y[first[i]:last[i]] = 0.0
      for (bl in itemList[[i]])
      {
        ii = which(booklet[[bl]]$first==first[i])
        g = mean_ElSym(b, a, booklet[[bl]]$first[-ii], booklet[[bl]]$last[-ii])
        lg1 = length(g) - 1
        scale_g=choose(lg1, 0:lg1)
        for (j in first[i]:last[i])
        {
          if (a[j] == 0) y[j]=y[j]+z[bl]*sum(g*scale_g*head(booklet[[bl]]$lambda,-a[last[i]]))
          if ((a[j] != a[last[i]])&(a[j]!=0)) y[j]=y[j]+z[bl]*sum(g*scale_g*head(tail(booklet[[bl]]$lambda,-a[j]),-(a[last[i]]-a[j])))
          if (a[j] == a[last[i]]) y[j]=y[j]+z[bl]*sum(g*scale_g*tail(booklet[[bl]]$lambda,-a[j]))
        }
      }
      b[first[i]:last[i]] = rgamma(1+last[i]-first[i],shape=sufI[first[i]:last[i]]+1.1,rate=y[first[i]:last[i]])
    }
    # identify 
      # within items
    for (i in 1:n)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
      # Between items
    f=b[2]
    b[-first] = b[-first]/f
    b[is.nan(b)] = 1 # deal with items that are not in data
      # Lambda
    for (bl in 1:nb) 
    {
      booklet[[bl]]$lambda = booklet[[bl]]$lambda*f^(0:sum(a[booklet[[bl]]$last]))
      booklet[[bl]]$lambda = booklet[[bl]]$lambda/booklet[[bl]]$lambda[1]
      lx[[bl]][iter,]=booklet[[bl]]$lambda
    }
    
    bx[iter,] = b
    setTxtProgressBar(pb, value=iter)
  }
  close(pb)
  OPCML_out=toOPLM(a,bx, first, last, H=NULL)
  return(list(a=a, b=bx,lambda=lx, beta.cml=OPCML_out$delta))
}

### Estimate lambda from CML
## arguments first, last and scoretab are for specific booklet
est_lambda <- function(b, a, first, last, scoretab)
{
  ms=length(scoretab)
  N=sum(scoretab)
  g=elsym(b,a,first,last)
  lx=vector("numeric",ms)
  lx=lx*NA  # deal with unobserved scores
  for (s in 1:ms)
  {
    if (scoretab[s]>0) {
      lx[s]=scoretab[s]/(N*g[s])
    }
  }
  return(lx)
}


#####################################################################
### Do CML
## Must be changed to handle the case where 0 categories does not occur
## Or category-scores (a) are not increasing
calibrateCML <- function(booklet, sufI, a, first, last, nIter) {
  nb = length(booklet)
  ni=length(first)
  b=rep(1,length(a))
  ic=b
  EsufI=sufI
  ref_cat=2
  
  ## 
  nn=0
  for (bl in 1:nb) nn=nn+booklet[[bl]]$m
  nn=nn*ni
  
  ## Implicit Equations  ###
  converged=FALSE
  iter=0
  pb = txtProgressBar(min=0, max=nIter)
  while ((!converged)&(iter<=nIter))
  {
    iter=iter+1
    EsufI=EsufI-EsufI
    for (bl in 1:nb)
    {
      EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
    }
    b=b*sufI/EsufI
    converged=(max(abs(sufI-EsufI))/nn<1e-04)
    setTxtProgressBar(pb, value=iter)
  }
  if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
  
  ## identification: 
  ## within items 
  for (i in 1:length(first))
  {
    range=first[i]:last[i]
    b[range]=b[range]/b[first[i]]
  }
  ## Between items; for now one item parameter set to zero(one)
  b[-first]=b[-first]/b[ref_cat]
  #const=log(b[2])
  #for (i in 1:length(first))
  #{
  #  range=(first[i]+1):last[i]
  #  b[range]=b[range]/exp(a[range]*const)
  #}
  
  
  ###  NR  ###
  H=matrix(0,length(a),length(a))
  convergence=FALSE
  while (!convergence)
  {
    iter=iter+1
    EsufI=EsufI-EsufI
    H=H-H
    for (bl in 1:nb)
    {
      EsufI = EsufI + E.STEP0(b,a,booklet[[bl]]$first+1,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
      H     = H     + H.STEP(b,a,booklet[[bl]]$first+1,booklet[[bl]]$last,booklet[[bl]]$scoretab)
    }
    # identify
    for (i in 1:length(first))
    {
      H[first[i],first[i]]=1
      EsufI[first[i]]=sufI[first[i]]
    }
    H[ref_cat,]=0
    H[,ref_cat]=0
    H[ref_cat,ref_cat]=1
    EsufI[ref_cat]=sufI[ref_cat]
    b=b*exp(solve(H,sufI-EsufI))
    convergence=(max(abs(EsufI-sufI))/nn<1e-10)
    setTxtProgressBar(pb, value=iter)
  }
  close(pb)
  
  lx=list()
  length(lx)=nb
  names(lx)=names(booklet)
  for (bl in 1:nb) lx[[bl]]=est_lambda(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
  
  OPCML_out=toOPLM(a,b, first, last, H=H)
  return(list(b=b,H=H,beta.cml=OPCML_out$delta,acov.cml=OPCML_out$cov_delta,lambda=lx))
}


### Change parameterization and normalization to produce OPCML output
toOPLM = function(a,b, first, last, H=NULL)
{
  ## for now remove zero category manually
  if (!is.null(H)) H=H[-first,-first]
  if (is.matrix(b)) b=b[,-first]
  if (is.vector(b)) b=b[-first]
  a=a[-first]
  new_first=first
  new_last=last-1
  for (i in 2:length(first))
  {
    ncat=last[i]-first[i]
    new_first[i]=new_last[i-1]+1
    new_last[i]=new_first[i]+ncat-1
  }
  first=new_first
  last=new_last
  logb=log(b)
  ########################
  
  ### Bayesian
  if (is.matrix(b))
  {
    k=ncol(b)
    delta=b
    for (r in 1:nrow(b))
    {
      for (i in 1:length(first))
      {
        delta[r,first[i]]=-logb[r,first[i]]/a[first[i]]
        if ((last[i]-first[i])>0)
        {
          tmp=(logb[r,(first[i]+1):last[i]]-logb[r,first[i]:(last[i]-1)])
          tmp=tmp/(a[first[i]:(last[i]-1)]-a[(first[i]+1):last[i]])
          delta[r,(first[i]+1):last[i]]=tmp
        }
      }
      delta[r,]=delta[r,]-mean(delta[r,]) ## mean centered
    }
    acov=cov(delta)
  }else{   ### CML; b is a single vector
    k=length(b)
    ## construct linear transformation
    DD=matrix(0,k,k)
    tel=1
    for (i in 1:length(first))
    {
      for (j in 1:(last[i]-first[i]+1))
      {
        if (j==1){
          DD[tel,tel]=-1/a[tel]
        }else
        {
          DD[tel,tel-1]=-1/(a[tel-1]-a[tel])
          DD[tel,tel]=1/(a[tel-1]-a[tel])
        }
        tel=tel+1
      }
    }
    CC=matrix(-1/k,k,k); diag(CC)=(k-1)/k
    AA=CC%*%DD
    ## calculate Deltaś and asymp. variance cov matrix
    delta=AA%*%logb
    if (!is.null(H))
    {
      acov=solve(H[-1,-1])
      acov=AA[,-1]%*%acov%*%t(AA[,-1])
    }else
    {
      acov=NULL
    }
  }  
  return(list(delta=delta, cov_delta=acov))
}

## Thus function expects categry thresholds delta, a vector of item_category scores a,
#  and first and last. All without the zero category.
#  It returns dexter parameters b, as well as new a, first and last with the zero category.
toDexter <- function(delta, a,first,last)
{
  ## Make D
  k=length(delta)
  DD=matrix(0,k,k)
  tel=1
  for (i in 1:length(first))
  {
    for (j in 1:(last[i]-first[i]+1))
    {
      if (j==1){
        DD[tel,tel]=-1/a[tel]
      }else
      {
        DD[tel,tel-1]=-1/(a[tel-1]-a[tel])
        DD[tel,tel]=1/(a[tel-1]-a[tel])        
      }
      tel=tel+1
    }
  }
  delta=delta-delta[1]  # normalize different
  b=exp(solve(DD)%*%delta) # exp(logb)
  names(b)=names(delta)
  
  new_first=first[1]
  new_last=last[1]+1
  for (i in 2:length(first))
  {
    nn=last[i-1]-first[i-1]
    new_first[i]=new_last[i-1]+1
    new_last=c(new_last,new_first[i]+nn+1)
  }
  new_a=vector("numeric",length(a)+length(first))
  new_b=vector("numeric",length(a)+length(first))
  new_a[-new_first]=a
  new_b[new_first]=1
  new_b[-new_first]=b
  
  ## put everything in a (minimal) parms object
  est=list(b=new_b, beta.cml=delta-mean(delta))
  inputs=list(ssIS=list(item_score=new_a),ssI=list(first=new_first,last=new_last))
  parms = list(est=est, inputs=inputs)
  return(parms)
}



############################ Elementary Symmetry
# All adapted for inclusion of zero-category
##
elsym <- function(b,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  g=double((ms+1))
  tmp = .C("ElSym",
          as.double(b),
          as.integer(a),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(-1),
          as.integer(-1),
          as.integer(n),
          as.integer(ms),
          as.double(g))
  return(tmp[[9]])
}


#################################### Means esf's
mean_ElSym <- function(b,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  g=double((ms+1))
  tmp = .C("meanElSym",
          as.double(b),
          as.integer(a),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(-1),
          as.integer(-1),
          as.integer(n),
          as.integer(ms),
          as.double(g),
          as.integer(0))
  return(tmp[[9]])
}



######################################## FUnctions for CML
H.STEP <- function(b,a,first,last,nscore)
{
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)^2)
  tmp=.C("H",
         as.double(b),
         as.integer(a),
         as.integer(length(b)),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[9]]
  tmp=as.matrix(tmp)
  dim(tmp)=c(length(b),length(b))
  tmp=tmp+t(tmp)
  diag(tmp)=diag(tmp)/2
  return(tmp)
}

#########################################
E.STEP0 <- function(b,a,first,last,nscore)
{
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)) 
  tmp=.C("E",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[8]]
  return(tmp)
}

########################################
E.STEP <- function(b,a,first,last,nscore)
{
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)) 
  tmp=.C("E_n",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[8]]
  return(tmp)
}



## produces a matrix of statistics for pairwise DIF
PairDIF_ <- function(par1,par2,cov1,cov2)
{
  labs=rownames(par1)
  D=kronecker(par2,t(par2),FUN="-")-kronecker(par1,t(par1),FUN="-") 
  var1=diag(cov1)
  var2=diag(cov2)
  S=(kronecker(var1,t(var1),FUN="+")-2*cov1)+(kronecker(var2,t(var2),FUN="+")-2*cov2)
  diag(S)=1
  D=D/sqrt(S)
  colnames(D)=labs; rownames(D)=labs
  return(D)
}

## produces a statistics for overall-DIF
OverallDIF_ <- function(par1,par2, cov1,cov2)
{
  r=1
  nI=length(par1)
  beta=par1-par2
  Sigma=cov1+cov2
  DIF_test=mahalanobis(beta[-r],rep(0,(nI-1)),Sigma[-r,-r])
  DIF_p=pchisq(DIF_test,(nI-1),lower.tail=FALSE)
  return(list(stat=DIF_test,df=nI-1, p=DIF_p))
}

###


#################################
bty = function (n, h = c(265, 75), c. = c(61, 66),
                l = c(25, 73), power = c(0.7, 1.742),
                fixup = TRUE, gamma = NULL, alpha = 1, ...)
{
  if (!is.null(gamma))
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L)
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- rep(c., length.out = 2L)
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, 0, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                       C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) *
                         rval), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}

##########################################
#' A print method for ENORM parms
#'
#' @param x An object produced by function \code{fit_enorm}
#' @param ... Any other parameters to the print method
#' @method print prms
#'
#'
print.prms <- function(x, ...){
  
  hpd=function(x, conf=0.95){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(paste0("(",as.character(round(x[ nnn ],digits=3))," , "
                  ,as.character(round(x[ n-nn+nnn ],digits=3)),")") )
  }
  
  if (x$inputs$method=="CML")
  {
    atab=as.data.frame(cbind(x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                             x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                             round(x$est$beta.cml,digits=3),
                             round(sqrt(diag(x$est$acov.cml)),digits=3)))
    colnames(atab)=c("item_id" ,"a", "B", "SE(B)")
  }

  if (x$inputs$method=="Bayes"){
    hh=apply(x$est$beta.cml,2,hpd)
    atab=as.data.frame(cbind(x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                             x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                             round(colMeans(x$est$beta.cml),digits=3),
                             round(apply(x$est$beta.cml, 2, sd),digits=3),
                             hh))
    colnames(atab)=c("item_id" ,"a", "mean(B)", "SD(B)", "95%hpd(B)")
  }
  row.names(atab)=NULL
  print(atab)
}

##########################################
#' Coerce parameters object to a data.frame of parameters
#'
#' @param x An object produced by function \code{fit_enorm}
#' @param ... Any other parameters to the as.data.frame
#' @method as.data.frame prms
#'
#'
as.data.frame.prms <- function(x, ...){
  hpd=function(x, conf=0.95){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(data.frame(l=round(x[ nnn ],digits=3),r=round(x[ n-nn+nnn ],digits=3)))
  }
  
  if (x$inputs$method=="CML")
  {
    atab=as.data.frame(cbind(x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                             x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                             round(x$est$beta.cml,digits=3),
                             round(sqrt(diag(x$est$acov.cml)),digits=3)))
    colnames(atab)=c("item_id" ,"a", "B", "SE(B)")
  }
  
  if (x$inputs$method=="Bayes"){
    hh=apply(x$est$beta.cml,2,hpd)
    atab=as.data.frame(item_id=x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                       a=x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                       mb=round(colMeans(x$est$beta.cml),digits=3),
                       sdb=round(apply(x$est$beta.cml, 2, sd),digits=3),
                       hpdl= hh[,1], hpdr=hh[,2])
    colnames(atab)=c("item_id" ,"a", "mean(B)", "SD(B)", "95%hpd(B) left","95%hpd(B) right")

  }
  row.names(atab)=NULL
  return(atab)
}




##########################################
#' A plot method for the interaction model
#'
#' Plot the item-total regressions fit by the interaction (or Rasch) model
#'
#'
#' @param x An object produced by function \code{fit_inter}
#' @param items The items to plot (column numbers). If NULL, all items will be plotted
#' @param summate If FALSE, regressions for polytomous items will be shown for each
#' response option separately; default is TRUE.
#' @param overlay If TRUE and more than one item is specified, there will be two plots,
#' one for the Rasch model and the other for the interaction model, with all items
#' overlayed; otherwise, multiple plots with the two models overlayed. Default is FALSE
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param show.observed If TRUE, the observed proportion correct at each sum score
#' will be shown as dots. Default is FALSE.
#' @param ... Any additional plotting parameters
#' @method plot rim
#'
plot.rim <- function(x, items=NULL, summate=TRUE, overlay=FALSE,
                    nc=1, nr=1, curtains=10, show.observed=FALSE, ...){
  allItems = x$ss$il$item_id
  qua = curtains/200
  if(qua>0 & qua<.5) {
    qnt = quantile(rep(as.integer(x$ss$tl$sumScore), x$ss$tl$N), c(qua,1-qua))
  } else {
    qnt=NULL
  }
  if (is.null(items)) items=allItems
  npic = length(items)
  if (length(items)==1) nr=nc=1
  if (overlay & !summate) overlay=FALSE
  if (overlay) {
  # only summate possible
  if (nr*nc==2) graphics::layout(matrix(1:2,nr,nc)) else graphics::layout(1)
  # do the Rasch model
  #
  z = x$regs$itrRM
  z = z[row.names(z) %in% items,]
  maxScore = ncol(z)-1
  graphics::plot(c(0,maxScore),c(0,max(z)),type="n",main="Rasch model",xlab="Test score",
                 ylab="Item score")
  if (!is.null(qnt)) {
    tmp = graphics::par('usr')
    graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
    graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
  }
  for (i in 1:npic) graphics::lines(0:maxScore,z[i,]) # the actual lines
  lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
  for (i in 1:npic) {
    graphics::points(lx[i], z[i,lx[i]+1], co="white", cex=1.6, pch=19)
    graphics::text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
  }
  # do the Interaction model
  #
  z = x$regs$itrIM
  z = z[row.names(z) %in% items,]
  maxScore = ncol(z)-1
  graphics::plot(c(0,maxScore),c(0,max(z,na.rm=TRUE)),type="n",main="Interaction model",xlab="Test score",
                 ylab="Item score")
  if (!is.null(qnt)) {
    tmp = graphics::par('usr')
    graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
    graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
  }
  for (i in 1:npic) graphics::lines(0:maxScore,z[i,]) # the actual lines
  lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
  for (i in 1:npic) {
    graphics::points(lx[i], z[i,lx[i]+1], co="white", cex=1.6, pch=19)
    graphics::text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
  }
  graphics::box()
  # end of overlay
  
  } else {
    # not overlay: do many plots
    ly = my_layout(npic, nr, nc)
    graphics::layout(matrix(1:(ly$nr*ly$nc), byrow=TRUE, ncol=ly$nc))
    if (summate) {
      # for each item in turn, do both models (with summation), and plot
      zI = x$regs$itrIM
      zR = x$regs$itrRM
      maxScore = ncol(zR)-1
      for (i in items) {
        mxY = max(zR[row.names(zR)==i,],na.rm=TRUE)
        graphics::plot(c(0,maxScore), c(0,mxY), type="n",
                       main=paste("Item", i),
                       xlab="Test score", ylab="Item score")
        if (!is.null(qnt)) {
          tmp = graphics::par('usr')
          graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
          graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
        }
        if(show.observed) {
          plt = x$ss$plt[x$ss$plt$item_id==i,]
          graphics::points(plt$sumScore,plt$meanScore,col="coral",pch=20)
        }
        graphics::lines(0:maxScore, zI[row.names(zI)==i,], col="gray80", lwd=3)
        graphics::lines(0:maxScore, zR[row.names(zR)==i,])
      }
      graphics::box()
    } else {
      zI = x$regs$ctrIM
      zR = x$regs$ctrRM
      maxScore = ncol(zR)-1
      # for each item in turn, similar but possibly multiline and coloured
      for (i in items) {
        prb = zI[x$ss$il$first[x$ss$il$item_id==i]:x$ss$il$last[x$ss$il$item_id==i],]
        pte = bty(nrow(prb), alpha=.6)
        graphics::plot(c(0,maxScore), 0:1, type="n",
                       main=paste("Item", i),
                       xlab="Test score", ylab="Probability")
        if (!is.null(qnt)) {
          tmp = graphics::par('usr')
          graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
          graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
        }
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j], lwd=3)
        }
        prb = zR[x$ss$il$first[x$ss$il$item_id==i]:x$ss$il$last[x$ss$il$item_id==i],]
        pte = bty(nrow(prb))
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j])
        }
      } # eol items
      graphics::box()
    } # eo not summate
  } # eo not overlay
}






##############################
my_layout <- function(npic, nr, nc) {
  if(npic==1) nr=nc=1
  nc = min(nc, 3)
  nc = min(nc, npic)
  nw = npic %/% nc + npic %% nc
  nr = min(nr, 3)
  nr = min(nr, nw)
  list(nr=nr, nc=nc)
}


##################################################
shinierInput <- function(FUN, id, namez, ...){
  inputs = character(length(namez))
  for (i in 1:length(namez))
  {
    inputs[i] = as.character(FUN(paste0(id, namez[i]), ...))
  }
  inputs
}



##################################
#' Derive scoring rules from keys
#'
#' For multiple choice items that will be scored as 0/1, derive the
#' scoring rules from the keys to the correct responses
#'
#'
#' @param keys  A data frame containing columns \code{item_id}, \code{nOptions}, and
#' \code{key} (the spelling is important). See details.
#' @return A data frame that can be used as input to \code{start_new_project}
#' @details
#' This function might be useful in setting up the scoring rules when all items
#' are multiple-choice and scored as 0/1. (Hint: Because the order in which the
#' scoring rules is not important, one can use the function to generate rules for
#' many MC items and then append per hand the rules for a few complex items.)
#'
#' The input data frame must contain the exact name of each item, the number
#' of options, and the key. If the keys are all integers, it will be assumed that
#' responses are coded as 1 through nOptions. If they are all uppercase letters,
#' it is assumed that responses are coded as A,B,C,... All other cases result
#' in an error.
#'
keys_to_rules <- function(keys) {
  # for backward compatibility we rename
  names(keys)[names(keys)=='item'] = 'item_id'
  
  if (is.numeric(keys$key)) ABC=FALSE else {
    if (all(keys$key %in% LETTERS)) ABC=TRUE
    else stop("You have inadmissible keys")
  }
  if (ABC) {
    m = match(keys$key, LETTERS)
    if (any(m>keys$nOptions)) stop("You have out-of-range keys")
    r = keys %>% group_by(.data$item_id) %>% do({
      y = data.frame(response=LETTERS[1:.$nOptions[1]], score=0L)
      y$score[match(.$key[1],LETTERS)] = 1
      y
    })
  } else {
    if (any(keys$key>keys$nOptions)) stop("You have out-of-range keys")
    r = keys %>% group_by(.data$item_id) %>% do({
      y = data.frame(response=1:.$nOptions[1], score=0L)
      y$score[.$key[1]] = 1
      y
    })
  }
  r
}

###########################################
#' A print method for the interaction model
#'
#' Print the available items for plots of the Rasch and the interaction models
#'
#'
#' @param x An object produced by function \code{fit_inter}
#' @param ... Included to stop check from nagging
#' @method print rim
#'
print.rim <- function(x, ...){
  available_items = x$ss$il[,"item_id"]
  print.default(available_items)
  invisible(x)
}

# use for forwarding arguments to e.g. plot function
merge_arglists = function(args, default = NULL, override = NULL)
{
  if(!is.null(default))
    for(nm in names(default)) 
      if(! nm %in% names(args)) 
        args[[nm]] = default[[nm]]
      
      if(!is.null(override))
        for(nm in names(override)) 
          args[[nm]] = override[[nm]]
        
        return(args)
}


#########################
#' Verbal aggression data
#' 
#' A data set of self-reported verbal behaviour in different frustrating
#' situations (Vansteelandt, 2000)
#' 
#' 
#' @name verbAggrData
#' @docType data
#' @format A data set with 316 rows and 25 columns.
#' @keywords datasets
NULL

###############################################
#' Scoring rules for the verbal aggression data
#' 
#' A set of (trivial) scoring rules for the verbal 
#' aggression data set
#' 
#' 
#' @name verbAggrRules
#' @docType data
#' @format A data set with 72 rows and 3 columns (item, response, score).
#' @keywords datasets
NULL

################################################
#' Item properties in the verbal aggression data
#' 
#' A data set of item properties related to the verbal
#' aggression data
#' 
#' 
#' @name verbAggrProperties
#' @docType data
#' @format A data set with 24 rows and 5 columns.
#' @keywords datasets
NULL

###########################################
#' Item properties in the PISA 2012 example
#' 
#' A data set of item properties in the PISA 2012 example (see the
#' help screen for function add_booklet)
#' 
#' 
#' @name PISA_item_class
#' @docType data
#' @format A data set with 109 rows and 6 columns.
#' @keywords datasets
NULL



