##################### Function relating to ability estimation
## Plausible Values
# pv1 <- function(b,a,first,last,score,mu=0,sigma=2)
# {
#   tmp=.C("PV1",
#          as.double(b),theta_
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

# Wrapper for C function that generates plausible values using a straightforward rejection algorithm.
# This one does not use recycling
# Normal prior(s)
# @param b          vector of beta's per item_score, including 0 category, ordered by item_id, item_score
# @param a          vector of discriminations per item_score, inclusing 0 category, ordered by item_id, item_score
# @param score      vector of sum scores
# @param npv        nbr of plausible values per sum score
# @param pop        vector of indices of population for each sum score
# @param mu, sigma  prior mu and sigma for each population
# @return
# matrix with in each row npv plausible values for each sum score
pv_ <- function(b,a,first,last,score,npv,mu,sigma,pop)
{
  tmp=.C("PV",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(mu),
         as.double(sigma),
         as.integer(score),
         as.integer(pop-1),
         as.integer(length(score)),
         as.integer(length(first)),
         as.integer(npv),
         as.double(rep(0*score,npv)))[[12]]
  tmp=as.vector(tmp)
  if (npv>1) dim(tmp)=c(length(score),npv)
  return(tmp)
}

# score is a vector of scores
# mu and sigma are scalar. One prior assumed
# returned is a matrix with for each score (=row) npv plausible values
pv_recycle <- function(b,a,first,last,score,npv,mu,sigma)
{
  # could be made slightly faster if score already sorted
  #stopifnot(npv==1)
  #pop = unique(pop)
  #stopifnot(length(pop)==1)
  
  if (sigma<0) stop('prior standard deviation must be positive')
  mx = sum(a[last])
  
  scrs = tibble(score=score) %>% 
    mutate(indx = row_number()) %>%
    arrange(.data$score)
  
  scr_tab = scrs %>%
    group_by(.data$score) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    right_join(tibble(score=0L:mx),by='score') %>%
    mutate(n=coalesce(n,0L)) %>%
    arrange(.data$score) %>%
    pull(.data$n) * npv
  
  cscr_tab = cumsum(scr_tab)

  theta = rep(0,length(score)*npv)
  

  tmp=.C("PVrecycle",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(mu),
         as.double(sigma),
         as.integer(scr_tab),
         as.integer(cscr_tab),
         as.integer(length(score)*npv),
         as.integer(length(first)),
         as.integer(max(a)),
         as.double(theta))[[12]]
  
  if(npv>1)
  {
    tmp = matrix(tmp, length(score), npv,byrow=TRUE)[order(scrs$indx),]
  } else
  {
    tmp = tmp[order(scrs$indx)]
  }
  
  return(tmp)
}


# Wrapper for C function that generates plausible values using a straightforward rejection algorithm
# Mixture of two normals as flexible prior
# @param b          vector of beta's per item_score, including 0 category, ordered by item_id, item_score
# @param a          vector of discriminations per item_score, inclusing 0 category, ordered by item_id, item_score
# @param score      vector of sum scores
# @param npv        nbr of plausible values per sum score
# @param p          prior group membership probabilities
# @param mu, sigma  mu and sigma for each population
# @param pop        vector of indices of population for each sum score. Currently not used
# @return
# matrix with in each row npv plausible values for each sum score
pv_mix <- function(b,a,first,last,score,npv,p,mu,sigma, pop)
{
  tmp=.C("PVMix",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(p),
         as.double(mu),
         as.double(sigma),
         as.integer(score),
         as.integer(length(score)),
         as.integer(length(first)),
         as.integer(npv),
         as.double(rep(0*score,npv)))[[12]]
  tmp=as.vector(tmp)
  if (npv>1) dim(tmp)=c(length(score),npv)
  return(tmp)
}

## Uses recycling to produce npv plausible values per score in scores
#  on a test defined by first and last. Scores defined to be 
#  0 to the maxScore on the test
#  If a prior sample is provided this is used; Otherwise prior is normal(mu,sigma)
#  A is a vector like a with alternative weightes. Used in theta_tables.
#  nscore is optional: a vector containing for each possible test score a frequency
#  is !is.null(nscore) and !is.null(A) I assume that nscore refers to scores with A as weights
#TO DO: This function needs to go in favour of pv_recycle
recycle_pv = function(b, a, first, last, npv=1, mu=0, sigma=2, nscore = NULL, prior_sample=NULL, A=NULL)
{
  if (is.null(A)){
    scores=0:sum(a[last])
    not_possible=setdiff(scores,possible_scores(a,first,last))
  }else {
    if (length(A)!=length(a)) stop("Wrong vector A in recycle_pv")
    if (any(A)<0) stop("Negative weights in vector A in recycle_pv")
    scores=0:sum(A[last])
    not_possible=setdiff(scores,possible_scores(A,first,last))
  }
  
  if (is.null(nscore)) nscore=rep(1,length(scores))
  n=nscore*npv
  if (length(not_possible)>0) n[not_possible+1]=0
  R=rep(0,sum(n))
  
  nI=length(first)
  if (is.null(A))
  {
    if (is.null(prior_sample))
    {
      if (sigma<0) stop('prior standard deviation must be positive')
      tmp=.C("recyclePV",
             as.double(b),
             as.integer(a),
             as.integer(first-1),
             as.integer(last-1),
             as.integer(nI),
             as.integer(n),
             as.double(mu),
             as.double(sigma),
             as.double(R))[[9]]
    }else
    {
      nprior=length(prior_sample)
      if (nprior<1) stop("wrong prior sample")
      tmp=.C("recyclePV2",
             as.double(b),
             as.integer(a),
             as.integer(first-1),
             as.integer(last-1),
             as.integer(nI),
             as.integer(n),
             as.double(prior_sample),
             as.double(nprior),
             as.double(R))[[9]]
    }
  }else
  {
    if (sigma<0) stop('prior standard deviation must be positive')
    tmp=.C("recyclePVaA",
           as.double(b),
           as.integer(a),
           as.integer(A),
           as.integer(first-1),
           as.integer(last-1),
           as.integer(nI),
           as.integer(n),
           as.double(mu),
           as.double(sigma),
           as.double(R))[[10]]
  }
  
  
    dim(tmp)=c(npv,length(scores))
    tmp=t(tmp)
    tmp[not_possible+1,]=NA
  
  return(tmp)
}

####################### FUNCTIONS TO UPDATE PRIOR of PLAUSIBLE VALUES

# Given samples of plausible values from one or more normal distributions
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values
#
# @param pv     vector of plausible values
# @param pop    vector of population indexes. These indices refer to mu and sigma
# @param mu     current means of each population
# @param sigma  current standard deviation of each population
#
# @return       a sample of mu and sigma

#TO DO: Adapt when there are multiple groups. Use hyperprior to keep means together
update_pv_prior<-function(pv, pop, mu, sigma)
{
  min_n=5 # no need to update variance when groups are very small
  for (j in 1:length(unique(pop)))
  {
    pv_group=subset(pv,pop==j)
    m=length(pv_group)
    if (m>min_n)
    {
      pv_var=var(pv_group)         
      sigma[j] = sqrt(1/rgamma(1,shape=(m-1)/2,rate=((m-1)/2)*pv_var))
    }
    pv_mean=mean(pv_group)
    mu[j] = rnorm(1,pv_mean,sigma[j]/sqrt(m)) 
  }
  return(list(mu=mu, sigma=sigma))
}


# Given samples of plausible values from a mixture of two normal distributions
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values. 
#
# Loosely based on Chapter 6 of 
# Marin, J-M and Robert, Ch, P. (2014). Bayesian essentials with R. 2nd Edition. Springer: New-York
# but works surpringly better than the official implementation in bayess
#
# @param pv     vector of plausible values
# @param p  vector of group membership probabilities
# @param mu     current means of each group
# @param sigma  current standard deviation of each group
# @param pop    vector of population indexes. Currently not used
#
# @return       a sample of p, mu and sigma
update_pv_prior_mixnorm = function (pv, p, mu, sigma, pop) {
  n = length(pv)
  z = rep(0, n)
  nj = c(0,0); 
  l = rep(1,2)
  v = rep(5,2)
  # Burnin and nr of Iterations with 0<Bin<nIter
  Bin = 1 
  nIter = 2

  mug = sig = pgr = matrix(0,nIter,2)
  mug[1,] = mu
  sig[1,] = sigma
  pgr[1,] = p
  mean_pv = mean(pv)
  var_pv = var(pv)
  for (i in 1:nIter)
  {
      ## latent group membership
    for (t in 1:n)
    {
      prob = pgr[i,]*dnorm(pv[t], mean = mug[i,], sd = sig[i,])
      if (sum(prob)==0) prob=c(0.5,0.5)
      z[t] = sample(1:2, size = 1, prob=prob)
    }
      ## Means
    for (j in 1:2)
    {
      nj[j]  = sum(z==j) 
      mug[i,j]  = rnorm(1, mean = (l[j]*mug[i,j]+nj[j]*mean(pv[z==j]))/(l[j]+nj[j]), 
                             sd = sqrt(sig[i,j]^2/(nj[j]+l[j])))
    }
      ## Vars 4=var_pv
    for (j in 1:2) 
    {
      if (nj[j]>1)
      {
        var_pvj = var(pv[z==j])
        sig[i,j] = sqrt(1/rgamma(1, shape = 0.5*(v[j]+nj[j]) ,
                                     rate = 0.5*(var_pv + nj[j]*var_pvj + 
                                           (l[j]*nj[j]/(l[j]+nj[j]))*(mean_pv - mean(pv[z==j]))^2)))
      }
    }
      ## membership probabilities
    p = rgamma(2, shape = nj + 1, scale = 1)
    pgr[i,] = p/sum(p)
  }
  return(list(p = colMeans(pgr[Bin:nIter,]), mu = colMeans(mug[Bin:nIter,]), sigma = colMeans(sig[Bin:nIter,]), grp = z))
}


# Plausible values. This one uses recycling 
# @param x                tibble(booklet_id <char or int>, person_id <char>, sumScore <int>, pop <int>)
# @param design           list: names: as.character(booklet_id), values: tibble(first <int>, last <int>) ordered by first
# @param b                vector of b's per item_score, including 0 category, ordered by item_id, item_score
# @param a                vector of weights per item_score, inclusing 0 category, ordered by item_id, item_score
# @param nPV              number of plausible values to return per person
### If the parameters are estimated with calibrate_Bayes:
# @param from             burn-in: first sample of b that is used
# @param by               every by-th sample of b is used. If by>1 auto-correlation is avoided.
# @param prior.dist       Prior distribution
#
# @return
# tibble(booklet_id <char or int>, person_id <char>, sumScore <int>, nPV nameless columns with plausible values)
pv = function(x, design, b, a, nPV, from = 20, by = 5, prior.dist = c("normal", "mixture"))
{
  prior.dist = match.arg(prior.dist)
  nPop = length(unique(x$pop))
  n_prior_updates = 10
  
  x$booklet_id = as.character(x$booklet_id)
  
  #if (prior.dist == "mixture") priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2), grp=sample(1:2, length(x$pop), replace=T, prob=c(0.5,0.5)))
  if (prior.dist == "mixture")
  {
    priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2))
    x$grp = sample(1:2, nrow(x), replace=T, prob=c(0.5,0.5))
  } else if (prior.dist == "normal")
  {
    priors = list(mu=rep(0,nPop), sigma=rep(4,nPop))
  }
  
  
  if (is.matrix(b))
  {
    which.pv = seq(from,(from-by)*(from>by)+by*nPV,by=by)
    nIter=max(which.pv)
    if (nrow(b)<nIter) stop(paste("at least", as.character(nIter), "samples of item parameters needed in function pv"))
    out_pv=matrix(0,length(x$sumScore),nPV)
 
    apv=1
    pb = txtProgressBar(min=0, max=nIter)
    for(iter in 1:nIter)
    {
      if (prior.dist == "mixture")
      {
        x = x %>% 
             group_by(.data$booklet_id, .data$grp) %>%
             mutate(PVX = pv_recycle(b[iter,], a, 
                                     design[[.data$booklet_id[1]]]$first, 
                                     design[[.data$booklet_id[1]]]$last, 
                                     .data$sumScore, 1, 
                                     priors$mu[.data$grp[1]], priors$sigma[.data$grp[1]])) %>%
             ungroup() 
        
        priors = update_pv_prior_mixnorm(x$PVX, priors$p, priors$mu, priors$sigma)
        x$grp = priors$grp
        
      } else if (prior.dist == "normal") 
      {
        x = x %>% 
          group_by(.data$booklet_id,.data$pop) %>%
          mutate(PVX = pv_recycle(b[iter,], a, 
                                  design[[.data$booklet_id[1]]]$first, 
                                  design[[.data$booklet_id[1]]]$last, 
                                  .data$sumScore, 1, 
                                  priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
          ungroup()
        
        priors = update_pv_prior(x$PVX, x$pop, priors$mu, priors$sigma)
      }

      if (iter == which.pv[apv])
      {
        colnames(x)[colnames(x)=='PVX'] = paste0('PV', iter)
        apv=apv+1
      }
      setTxtProgressBar(pb, value=iter)
    }
    close(pb)
    return( select(x, .data$booklet_id, .data$person_id,.data$sumScore, matches('PV\\d+')))
    
  }else
  {
   for(iter in 1:n_prior_updates) 
   {

     pv = x %>%
       group_by(.data$booklet_id,.data$pop) %>%
       mutate(pv = pv_recycle(b, a, 
                              design[[.data$booklet_id[1]]]$first, 
                              design[[.data$booklet_id[1]]]$last, 
                              .data$sumScore, 1, 
                              priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
       ungroup() 

      priors = update_pv_prior(pv$pv,pv$pop,priors$mu,priors$sigma)
   }

  return(  
     x %>% 
       group_by(.data$booklet_id, .data$pop) %>%
       do({
         bkID = .$booklet_id[1]
         popnbr = .data$pop[1]
         out_pv = pv_recycle(b, a, design[[bkID]]$first, design[[bkID]]$last, .$sumScore, nPV, priors$mu[popnbr], priors$sigma[popnbr])
         data.frame(.$person_id, .$sumScore, as.data.frame(out_pv), stringsAsFactors = FALSE)
        }) %>%
       ungroup() %>%
       select(-.data$pop)) 
  }
}

# simulate 0-based index of responses to a single item. Adapted for inclusion zero
# the index i is local in that it refers to an item-index in first and last
renorm <- function(b,a,theta,first,last,i)
{
  m=length(theta)
  x=rep(0,m)
  if ((i<0)|(i>length(first))) stop("wrong item-indx in renorm")
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

# simulate responses to a single item. 
# NOt adapted for inclusion of parameter for 
# zero category
renorm0 <- function(b,a,theta,first,last,i)
{
  m=length(theta)
  x=rep(0,m)
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

#void sampleNRM2(double *theta, double *b, int *a, int *first, int *last, int *nI, int *m, int *test_score)
# simulate test-scores rather then response patterns. Adapted for inclusion zero
rscore <- function(theta,b,a,first,last, cntr=NULL, use_b_matrix=FALSE)
{
  if(use_b_matrix)
  {
    if(is.null(cntr)) stop('use_b_matrix is true, need a counter')
    b = b[cntr$get(),]
  }
  m=length(theta)
  nI=length(first)
  x=rep(0,m)
  tmp=.C("sampleNRM2",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nI),
         as.integer(m),
         as.integer(x))[[8]]
  return(tmp)
}


# Expected scores given a single ability value theta
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



# Expected distribution given one ability theta
pscore <- function(theta, b, a, first, last)
{
  g=elsym(b, a, first, last)
  p=rep(0, length(g))
  for (s in 1:length(g))
  {
    p[s]=g[s]*exp((s-1)*theta)
  }
  return(p/sum(p))
}

####################################################
# Computes likelihood and test information for internal use
#
# For a vector of thetas it returns:
# l = a matrix (n of response cats * length of theta) of the likelihood or log-likelihood if log=TRUE
# I = a vector of the information function computed at each theta = sum(P'^2/PQ)
# J = something sum(P'P"/PQ) 
# The vector theta can be a set of quadrature points or 
# estimates to compute their SE
#
# Note: can not deal with Inf or NA values in theta
IJ_ = function(b,a,first, last, theta, log=FALSE)
{
  nI = length(first)
  nT = length(theta)
  max_ncat = max(last-first) + 1L
  I = rep(0, nI * nT)
  J = rep(0, nI * nT)
  logFi = rep(0, nT)
  
  tmp=.C("IJ_c",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(I),
         as.double(J),
         as.double(logFi),
         as.integer(nI),
         as.integer(nT),
         as.integer(max_ncat))

  scores = 0:sum(a[last])
  l = sweep(outer(scores,theta), 2, tmp[[8]], '-')
  if (!log) l=exp(l)
  return(list(I=colSums(matrix(tmp[[6]],nI,nT)), J=colSums(matrix(tmp[[7]],nI,nT)), l=l))
}


#### ML estimation of theta
# uses an implicit equations algorithm
# se is the option to calculate standard-errors or not
theta_MLE <- function(b,a,first,last, se=FALSE)
{
  mx=sum(a[last])
  theta=rep(0,mx-1)
  max_a = max(a)

  n=length(first)
  theta = .C("theta_mle_c",
           as.double(theta),
           as.double(b),
           as.integer(a),
           as.integer(first-1L),
           as.integer(last-1L),
           as.integer(n),
           as.integer(mx),
           as.integer(max_a))[[1]]

  sem = NULL
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem = c(NA, 1/sqrt(f$I), NA)
  }
  theta=c(-Inf,theta,Inf)
  
  return(list(theta=theta,se=sem))
}



##### EAP based on npv (default 500) plausible values ####
# Uses recycling to get npv plausible values for each sum score.
# Options:
# Smoothing is an option to avoid non-monotonicity due to random fluctuations
# se is the option to calculate standard-errors or not
### In R-code:
# n=rep(npv,length(score))
# R=matrix(0,length(score),npv)
# while (any(n>0))
# {
#   atheta=rnorm(1,mu,sigma)
#   sm=rscore(b,a,atheta,first,last)+1
#   if (n[sm]>0)
#   {
#     R[sm,n[sm]]=atheta
#     n[sm]=n[sm]-1
#   }
# }
# The final argument allows EAPs based on A-weighted score, where A need not equal a.
####
theta_EAP <- function(b, a, first, last, npv=500, mu=0, sigma=4, smooth=FALSE, se=FALSE, A=NULL)
{
  if (!is.null(A)){
    tmp=recycle_pv(b, a, first, last, npv=npv, mu=mu, sigma=sigma, A)
    mx = sum(A[last])
    theta = rep(NA,(mx+1))
  }else
  {
    score = possible_scores(a,first,last)
    tmp = pv_recycle(b,a,first,last,score,npv,mu,sigma)
    mx = sum(a[last])
    theta = rep(NA,(mx+1))
  }
  theta[score+1]=rowMeans(tmp)
  if (is.null(A))
  {
    if (smooth) {
      score=0:mx
      theta = predict(lm(theta ~ poly(score,7)))
    }
    sem=rep(NA,(mx+1))
  }else
  {
    if (smooth) {
      score=0:mx
      theta = predict(lm(theta ~ poly(score,7)))
    }
    sem=rep(NA,(mx+1))
  }
  if (se) sem=apply(tmp,1,sd)
  return(list(theta=theta, se=sem))
}

## EAP using Jeffrey's prior: aka propto sqrt(information)
# Uses a weighted average to integrate over a grid defined by:
# grid_from, grid_to and grid_length.
theta_jEAP=function(b, a, first,last, se=FALSE, grid_from=-6, grid_to=6, grid_length=101) 
{
  theta_grid = seq(grid_from, grid_to, length=grid_length)
  f = IJ_(b,a,first,last,theta_grid)
  prior=sqrt(f$I)
  w = sweep(f$l, 2, prior, '*')
  theta = apply(w, 1, function(x)weighted.mean(theta_grid, w=x))
  sem=rep(NA,length(theta))
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  return(list(theta=theta,se=sem))
}

# computes distribution of score weighted with A conditionally on score weighted with a
# used for estimation of theta from the unweighted score (for instance)
# colSums(G[,,(1-idx)+1]) are elementary symmetric function with a
elsymat <- function(b, a, A, first, last)
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
  #return(c(-Inf,theta,Inf))
  theta = c(-Inf,theta,Inf)
  return(list(theta=theta,se=rep(NA, length(theta))))
}


# Estimate a single ability for a whole score distribution.
# Used for the 3DC standard setting method
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


########################################################
## Score-by-score table. Currently using mean_ElSym
########################################################
# @param m        a rim object produced by fit_inter (but not yet documented anywhere what that looks like)
# @param AB       list: two mutually exclusive subsets of items as indexes of m$ss$il
# @param model    Character: Indicates which model is used: "Rasch" or "IM"
# @return         A list with tbl being a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b
# @details        NA's indicate that a total scores was not possible given the weights
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
  nI=length(first)
  bA=NULL; bB=NULL
  aA=NULL; aB=NULL
  cA=NULL; cB=NULL
  lastA=NULL; firstA=NULL
  lastB=NULL; firstB=NULL
  telAF=1; telBF=1
  for (i in 1:nI)
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
      if (g[s+1]>0)
      {
        for (s_a in max(0,s-MscB):min(s,MscA))
        {
          s_b=s-s_a
          out[s_a+1,s_b+1] = -Inf
          if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
          {
            out[s_a+1,s_b+1] = log(gA[s_a+1]) + log(gB[s_b+1]) - log(g[s+1])
            out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
          }
        }
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
      if (g[s+1]>0)
      {
        for (s_a in max(0,s-MscB):min(s,MscA))
        {
          s_b=s-s_a
          out[s_a+1,s_b+1] = -Inf
          if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
          {
            out[s_a+1,s_b+1] = log(gA[s_a+1])+log(gB[s_b+1])-log(g[s+1])
            out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
          }
        }
      }
    }
  }
  return(list(tbl=exp(out),m=m,AB=AB,model=model))
}

########################################################
## Score-by-score table ENORM Currently using mean_ElSym
########################################################
# @param m        a parms object produced by fit_enorm
# @param AB       list: two mutually exclusive subsets of items as indexes of items in first and last
# @return         A list with tbl being a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b. 
# @details        NA's indicate that a total scores was not possible given the weights
SSTable_ENORM <- function(b, a, first, last, AB) {
  A = AB[[1]]
  B = AB[[2]]
  ### Check
  if (length(intersect(A,B))!=0) stop("sets not disjunct")
  unionAB = unlist(AB, use.names = F)
  firstA =first[A]; lastA = last[A]
  firstB =first[B]; lastB = last[B]
  MscA = sum(a[lastA])
  MscB = sum(a[lastB])
  Msc  = sum(a[last[unionAB]])
  out = matrix(NA,MscA+1,MscB+1)
  
  gA = mean_ElSym(b,a,firstA,lastA)
  gB = mean_ElSym(b,a,firstB,lastB)
  g =  mean_ElSym(b,a,first[unionAB],last[unionAB])
  
  for (s in 0:Msc)
  {
    if (g[s+1]>0)
    {
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b = s-s_a
        out[s_a+1,s_b+1] = -Inf
        if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
        {
          out[s_a+1,s_b+1] = log(gA[s_a+1]) + log(gB[s_b+1]) - log(g[s+1])
          out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
        }
      }
    }
  }
  return(list(tbl=exp(out),AB=AB,model="RM"))
}


#################################### Item-Total Regressions


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

## low-level check for trivial scores. 
# If TRUE it means that there is no data-reduction
# Each weighted score can be obtained with only one response pattern
# @pi_mat An item-total regression matrix as produced by ittotmat
check_trivial_scores_<-function(pi_mat,scoretab=NULL)
{
  if (!is.null(scoretab)) pi_mat=pi_mat[,scoretab>0]
  pi_mat[(pi_mat==1)|(pi_mat==0)]=NA
  return(all(is.na(pi_mat)))
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


######################################################################
# Estimation of Rasch and Interaction Model 
#######################################################################
# Using Newton-Raphson per item
# Currently with elsym and scale
# @param ss list containing:
#   il: one row per item, ordered item_id
#       tibble(first, last, sufC <sum(item_score*sumScore)>, nCat <nbr of score categories including 0>)
#   sl: one row per item-scorecat(including 0), ordered item_id, item_score
#       tibble(item_score, sufI <count for item_score>, sufC <sum(item_score * sumScore)>)
#   tl: one row per testscore (only one test allowed), complete range 0:max_observed, ordered by test_score
#       tibble(N <count for test score>)
# @returns:
#      bRM:    Parameter estimates of the Rasch model
#      bIM:    Parameter estimates of the Interaction model
#      cIM:    Estimate of (log-)interaction parameter
#      cRM:    Interaction parameters under the Rasch model: all equal to 1
#      HRM:    Block diag. Asymptotic var-covvar matrix of item parameters under RM
#     se.c:   Standard error of interaction parameter
# fit.stat: log(cIM)/se.c. Wald statistic normally distributed under Rasch model
#########################################################################
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
  ic=rep(1,nI)
  var.ic=vector("numeric", nI)
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
  
  if (check_trivial_scores_(ittotmat0(b,ic[C],a,first,last),scoretab)){
    warning("Only trivial weighted scores observed") }
  
  converged=2
  scale=2
  while(converged>0.01)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last) ##
    pi_mat[is.na(pi_mat)]=0
    for (i in 1:nI)
    {
      if (length(upd_set[[i]])>0)
      {
        pi=pi_mat[upd_set[[i]],,drop=FALSE]
        E=sufI[upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%tcrossprod(diag(scoretab), pi) #diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # NR update for parameters of item i
        update=solve(H*scale,E)
        b[upd_set[[i]]]=b[upd_set[[i]]]*exp(update)
        converged=max(converged,max(abs(E))/m) #
        HRM[upd_set[[i]],upd_set[[i]]]=H
      }
    }
    if (converged<1) scale=1
  }

  bRM=b
  cRM=ic
  
  ## IM
  converged=2
  scale=2
  while(converged>0.001)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last)
    pi_mat[is.na(pi_mat)]=0
    for (i in 1:nI)
    {
      # gradient and hessian for thresholds of item i
      if (length(upd_set[[i]])>0)
      {
        pi=pi_mat[upd_set[[i]],,drop=FALSE]
        E=sufI[upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%tcrossprod(diag(scoretab), pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # gradient and hessian for interaction parameter
        ncol_pi=ncol(pi); nrow_pi=nrow(pi)
        E=c(E,sufC[i])
        H=cbind(H,rep.int(0,nrow(H)))
        H=rbind(H,rep.int(0,ncol(H)))
        k=1
        e0=0; e1=0
        f=matrix(0,nrow_pi,ncol_pi)
        g=matrix(0,nrow_pi,ncol_pi)
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
        var.ic[i]=solve(H)[nrow(H),nrow(H)]
        converged=max(converged,max(abs(E))/m)
      }
    }
    if (converged<1) scale=1
  }
  
  return(list(group=ss$group,bRM=bRM,cRM=cRM,bIM=b,cIM=ic,se.c=sqrt(var.ic),HRM=HRM, fit.stats=log(ic)/sqrt(var.ic)))
}



##################################################### calibrate incomplete designs: 
# CML
# Bayes
#####################################################

#### Bayes
calibrate_Bayes = function(itemList, booklet, sufI, b, a, first, last, nIter, fixed_b=NULL, lambda_out=FALSE) 
{
  nb = length(booklet)
  n = length(itemList)
  y = rep(0, length(sufI))
  z = NULL
  bx = matrix(0, nIter, length(b))
  lx=list()
  if (lambda_out)
  {
    length(lx)=nb
    names(lx)=names(booklet)
    for (bl in 1:nb) lx[[bl]]=matrix(0,nIter,sum(a[booklet[[bl]]$last])+1)
  }
 
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
      booklet[[bl]]$lambda[idx] = rgamma(length(idx), shape=booklet[[bl]]$scoretab[idx]+0.1, rate=(g*scale_g*z[bl])[idx]) # 1.1
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
      b[first[i]:last[i]] = rgamma(1+last[i]-first[i],shape=sufI[first[i]:last[i]]+1.1,rate=y[first[i]:last[i]]) #1.1
    }
    # identify
    for (i in 1:n)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    if (is.null(fixed_b))
    {
      f=b[2]
      b[-first] = b[-first]/f
      # Lambda
      for (bl in 1:nb) 
      {
        booklet[[bl]]$lambda = booklet[[bl]]$lambda*f^(0:sum(a[booklet[[bl]]$last]))
        booklet[[bl]]$lambda = booklet[[bl]]$lambda/booklet[[bl]]$lambda[1]
        if (lambda_out) lx[[bl]][iter,]=booklet[[bl]]$lambda
      }
    }else
    {
      fixed_set=which(!is.na(fixed_b))
      b[fixed_set]=fixed_b[fixed_set]
      if (lambda_out) { for (bl in 1:nb) lx[[bl]][iter,]=booklet[[bl]]$lambda }
    }
    b[is.nan(b)] = 1 # deal with items that are not in data
    bx[iter,] = b
    setTxtProgressBar(pb, value=iter)
  }
  close(pb)

  OPCML_out=toOPLM(a,bx, first, last, H=NULL,fixed_b=fixed_b)
  return(list(a=a, b=bx,lambda=lx, beta.cml=OPCML_out$delta))
}

### Estimate lambda from CML
## arguments first, last and scoretab are for specific booklet
est_lambda <- function(b, a, first, last, scoretab)
{
  ifelse(scoretab>0, scoretab/(elsym(b,a,first,last)*sum(scoretab)), NA) 
}  

  
#####################################################################
### Do CML
## Must be changed to handle the case where 0 categories does not occur
## Or category-scores (a) are not increasing
## If fixed_b is not NULL unfixed c.q. free parameters are NA
## Note that fixed_b contains values in the dexter parametrisation (including 0 category)

calibrate_CML <- function(booklet, sufI, a, first, last, nIter, fixed_b=NULL) {
  nb = length(booklet)
  ni = length(first)
  EsufI = sufI
  max_nr_iter = 30
  
  if (is.null(fixed_b)) # if no fixed parameters
  {
    nn= sum(sufI)
    b = rep(1,length(a))
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
      b = b*sufI/EsufI
      converged=(max(abs(sufI-EsufI))/nn<1e-04)
      setTxtProgressBar(pb, value=iter)
    }
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
        ### identification ###
    # within items
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    # between items. Note if ref_cat is allowed to be something else then 2 has to adapt toOPLM and toDexter
    ref_cat=2
    b[-first] = b[-first]/b[ref_cat]
    
              ###  NR  ###
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=3
    while ((!converged)&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI=EsufI-EsufI
      H=H-H
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
        H     = H     + H.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
      }
      # identify
      for (i in 1:ni)
      {
        H[first[i],first[i]]=1
        EsufI[first[i]]=sufI[first[i]]
      }
      H[ref_cat,]=0
      H[,ref_cat]=0
      H[ref_cat,ref_cat]=1
      EsufI[ref_cat]=sufI[ref_cat]
      b = b*exp(solve(H*scale,sufI-EsufI))
      converged=(max(abs(EsufI-sufI))/nn<1e-10)
      setTxtProgressBar(pb, value=iter)
      if (nr_iter==2) scale=1
    }
    close(pb)
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  }else  ### if fixed parameters
  {
    fixed_set=which(!is.na(fixed_b))
    update_set=which(is.na(fixed_b))
    b=fixed_b
    ni_free=sum(is.na(fixed_b[last]))
    b[update_set]=1
    nn=0
    for (bl in 1:nb) nn=nn+booklet[[bl]]$m
    nn=nn*ni_free
    
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
      b[update_set] = b[update_set]*sufI[update_set]/EsufI[update_set]
      converged=(max(abs(sufI[update_set]-EsufI[update_set]))/nn<1e-04)
      setTxtProgressBar(pb, value=iter)
    }
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=2
    while ((!converged)&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI=EsufI-EsufI
      H=H-H
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
        H     = H     + H.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
      }
      # identify
      for (i in 1:length(first))
      {
        H[first[i],first[i]]=1
        EsufI[first[i]]=sufI[first[i]]
      }
      H[fixed_set,]=0
      H[,fixed_set]=0
      diag(H)[fixed_set]=1
      EsufI[fixed_set]=sufI[fixed_set]
      b = b*exp(solve(H*scale,sufI-EsufI))
      converged=(max(abs(EsufI[update_set]-sufI[update_set]))/nn<1e-10)
      setTxtProgressBar(pb, value=iter)
      scale=1
    }
    close(pb)
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  }
  
  lx=list()
  length(lx)=nb
  names(lx)=names(booklet)
  for (bl in 1:nb) lx[[bl]]=est_lambda(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
  
  OPCML_out=toOPLM(a, b, first, last, H=H, fixed_b=fixed_b)
  return(list(b=b, H=H, beta.cml=OPCML_out$delta, acov.cml=OPCML_out$cov_delta, lambda=lx, n_iter=iter))
}

## Get the score distribution of a booklet from fit_enorm
#  If also gives you a smooth version based on a polynomial smoothing
#  of the log-lambda's
#### Example
# db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender=""))
# add_booklet(db, verbAggrData, "agg")
# f=fit_enorm(db)
# xx=dexter:::ENORM2ScoreDist(f,degree=7,booklet_id="data")
# plot(xx$score,xx$p.obs,pch=16, main=paste("degree =", as.character(xx$degree)))
# lines(xx$score,xx$p.obs,col="red")
# lines(xx$score,xx$p.smooth, col="green")
ENORM2ScoreDist <- function (parms, degree=7, booklet_id) 
{
  bk_indx=which(names(parms$inputs$bkList)==booklet_id)
  stopifnot(length(bk_indx)>0) 
  if (parms$inputs$method!="CML") stop("ENORM2ScoreDist only for CML at the moment")
  score_range=0:(length(parms$est$lambda[[bk_indx]])-1)
  a=parms$inputs$ssIS$item_score
  first=parms$inputs$bkList[[bk_indx]]$first
  last=parms$inputs$bkList[[bk_indx]]$last
  lambda=parms$est$lambda[[bk_indx]]
  b=parms$est$b
  log_l=log(lambda)
  
  degree=min(degree,length(!is.na(log_l)))
  qr=lm(log_l~poly(score_range,degree,raw=TRUE)) 
  beta=as.numeric(qr$coefficients)[-1]
  
  lambda[is.na(lambda)]=0
  mx = sum(a[last])
  g = elsym(b,a,first,last)
  sc_obs=vector("numeric", mx+1)
  sc_sm=vector("numeric", mx+1)
  num_obs=0.0
  num_sm=0.0
  for (s in 0:mx)
  {
    sc_obs[s+1]=g[s+1]*lambda[s+1]
    sc_sm[s+1]=g[s+1]*exp(sum(beta*s^(1:degree)))
    num_obs=num_obs+sc_obs[s+1]
    num_sm=num_sm+sc_sm[s+1]
  }
  sc_obs=sc_obs/num_obs
  sc_sm=sc_sm/num_sm
  data.frame(score=score_range,
             n.obs=sc_obs*sum(parms$inputs$stb$N), 
             n.smooth=sc_sm*sum(parms$inputs$stb$N),
             p.obs=sc_obs,
             p.smooth=sc_sm)
}


### Change parameterization and normalization to produce OPCML output
## assumes that the first parameter is the reference unless there are fixed parameters
toOPLM = function(a, b, first, last, H=NULL, fixed_b=NULL)
{
  ## for now remove zero category manually
  if (!is.null(H)) H=H[-first,-first]
  if (is.matrix(b)) {
    b=b[,-first]
    if (is.null(dim(b))) b=as.matrix(t(b))
  }
  if (is.vector(b)) b=b[-first]
  if (!is.null(fixed_b)) fixed_b=fixed_b[-first]
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
  acov=NULL
  ########################
  
  ### Bayesian # a better criterion may be needed if b is Bayesian but there was one sample and b is a vector
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
      if (is.null(fixed_b)) delta[r,]=delta[r,]-mean(delta[r,]) ## mean centered
    }
    if (nrow(delta)>20*ncol(delta)){
      acov=cov(delta)
    }
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
    if (is.null(fixed_b))
    {
      CC=matrix(-1/k,k,k); diag(CC)=(k-1)/k
      AA=CC%*%DD
    ## calculate Deltaś and asymp. variance cov matrix
    #  Note: assumes that the first parameters is the reference
      delta=AA%*%logb
      if (!is.null(H))
      {
        acov=solve(H[-1,-1])
        acov=AA[,-1]%*%acov%*%t(AA[,-1])
      }
    }else # if there are fixed parameters we do not normalize
    {
      delta=DD%*%logb
      if (!is.null(H))
      {
        fixed_set=which(!is.na(fixed_b))
        acov=solve(H[-fixed_set,-fixed_set])
        acov=DD[,-fixed_set]%*%acov%*%t(DD[,-fixed_set])
      }
    }
  }  
  return(list(delta=delta, cov_delta=acov, a=a, first=new_first, last=new_last))
}

## Thus function expects category thresholds delta, a vector of item_category scores a,
#  and first and last. All without the zero category.
#  It returns dexter parameters b, as well as new a, first and last with the zero category.
toDexter <- function(delta, a, first, last, re_normalize=TRUE)
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
  if (re_normalize) delta=delta-delta[1]  # normalize different
  b=exp(solve(DD)%*%delta) # exp(logb)
  names(b)=names(delta)
  
  new_first=first[1]
  new_last=last[1]+1
  if (length(first)>1)
  {
    for (i in 2:length(first))
    {
      nn=last[i]-first[i]
      new_first[i]=new_last[i-1]+1
      new_last=c(new_last,new_first[i]+nn+1)
    }
  }
  new_a=vector("numeric",length(a)+length(first))
  new_b=vector("numeric",length(a)+length(first))
  new_a[-new_first]=a
  new_b[new_first]=1
  new_b[-new_first]=b
  
  ## put everything in a (minimal) parms object
  est=list(b=new_b, a=new_a, beta.cml=delta-mean(delta))
  inputs=list(ssIS=list(item_score=new_a),ssI=list(first=new_first,last=new_last))
  parms = list(est=est, inputs=inputs)
  return(parms)
}

#################################### Wrappers for C- FUnctions for CML
## THis version uses Elsym and is adapted for inclusion of b_0. However
## THere must still be a small error somewhere.
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

## This version uses Elsym0... 
H.STEP0 = function(b,a,first,last,nscore)
{
  first=first+1
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)^2)
  tmp=.C("H0",
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

########################################
# version with Elsym
E.STEP <- function(b,a,first,last,nscore)
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
# version with Elsym0
E.STEP0 <- function(b,a,first,last,nscore)
{
  first=first
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)) 
  tmp=.C("E0",
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



