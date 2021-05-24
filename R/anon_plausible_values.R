#############################################################################
# R functions to generate plausible values
############################################################################


# score is a vector of scores
# if alpha > 0, mu and sigma should have length 2
# returned is a matrix with for each score (=row) npv plausible values
pv_recycle = function(b, a, first, last, score, npv, mu, sigma, alpha=-1, A=NULL)
{
  if (any(sigma<0)) stop('prior standard deviation must be positive')
  if(alpha<=0 && (length(mu)!=1 || length(sigma)!=1)) stop('alpha<0, mu and sigma should be length 1')
  if(alpha>0 && (length(mu)!=2 || length(sigma)!=2)) stop('alpha>0, mu and sigma should be length 2')
  
  a = as.integer(a)
  score = as.integer(score)
  
  scrs = tibble(score=score) %>% 
    mutate(indx = row_number()) %>%
    arrange(.data$score)
  

  scr_tab = as.integer(npv) * score_tab_single(score, sum(a[last]))


  if(is.null(A))
    A = a
  
  A = as.integer(A)
  
  # scr_tab is consumed in the process
  # so should return with all zeros
  theta = PVrecycle(b, a, as.integer(first-1L),as.integer(last-1L),
                          mu, sigma, scr_tab, A, alpha)[,1]

  if(npv>1)
  {
    theta = matrix(theta, length(score), npv, byrow=TRUE)[order(scrs$indx),,drop=FALSE]
  } else
  {
    theta = theta[order(scrs$indx)]
  }
  return(theta)
}

# when score is pre-sorted this function is equivalent to pv_recycle but slightly faster for large score vectors
# when score is not sorted, outputted plausible values will not match the order of the inputted scores
# however, all summary values will of course be correct so useful for prior updating
# if scr_tab is not null, some time is saved by using precomputed scr_tab
pv_recycle_sorted = function(b,a,first,last,score,npv,mu,sigma,alpha=-1,A=NULL, scr_tab=NULL)
{
  if (any(sigma<0)) stop('prior standard deviation must be positive')
  if(alpha<=0 && (length(mu)!=1 || length(sigma)!=1)) stop('alpha<0, mu and sigma should be length 1')
  if(alpha>0 && (length(mu)!=2 || length(sigma)!=2)) stop('alpha>0, mu and sigma should be length 2')
  
  if(length(score)==0)
    return(double())
  
  
  a = as.integer(a)

  if(is.null(scr_tab))
  {
    scr_tab = score_tab_single(score, sum(a[last]))
  }
  
  scr_tab = as.integer(npv) * scr_tab

  if(is.null(A))
    A = a
  
  A = as.integer(A)
  
  # scr_tab is consumed in the process
  # so should return with all zeros
  theta = PVrecycle(b, a, as.integer(first-1L),as.integer(last-1L),
                    mu, sigma, scr_tab, A, alpha)[,1]
  
  if(npv>1)
    return(matrix(theta, length(score), npv, byrow=TRUE))

  theta
}



####################### FUNCTIONS TO UPDATE PRIOR of PLAUSIBLE VALUES #####

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


update_pv_prior = function(pv, pop, mu, sigma)
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

# Update prior for plausible values as an hierarchical normal model
#
# @param pv       vector of plausible values
# @param pop      vector of group membership indices
# @param mu       vector of current group means
# @param sigma    vector of current standard deviations; assumed equal in each group
# @param mu.a     current overall mean of group means
# @param sigma.a  current standard deviation of group means
#
# @return       a sample of p, mu, sigma, mu.a, and sigma.a
#
# @details      Given samples of plausible values from an hierarchical model with >=2 groups
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values. Our implementation is based on code 
# provided by Gelman, A., & Hill, J. (2006). Data analysis using regression and 
# multilevel/hierarchical models. Cambridge university press.
#
# to do: test n>=5, anders geen sd update
# to do: Allow within group variance to be different
update_pv_prior_H = function(pv, pop, mu, sigma, mu.a, sigma.a)
{
  J =length(unique(pop))
  n = length(pv)
  sigma=sigma[1]
  for (j in 1:J){
    n.j = sum(pop==j)
    y.bar.j = mean(pv[pop==j])
    V.a.j = 1/(n.j/sigma^2 + 1/sigma.a^2)
    a.hat.j = V.a.j*((n.j/sigma^2)*y.bar.j + (1/sigma.a^2)*mu.a)
    mu[j] = rnorm(1, a.hat.j, sqrt(V.a.j))
  }
  mu.a = rnorm (1, mean(mu), sigma.a/sqrt(J))
  sigma = sqrt(sum((pv-mu[pop])^2)/rchisq(1,n-1))
  sigma.a = sqrt(sum((mu-mu.a)^2)/rchisq(1,J-1))
  
  return(list(mu=mu, sigma=rep(sigma,J), mu.a=mu.a, sigma.a=sigma.a))
}


# Update prior for plausible values as a mixture of two normals
#
# @param pv     vector of plausible values
# @param p      vector of group membership probabilities
# @param mu     current means of each group
# @param sigma  current standard deviation of each group
# @param prior.dist  Hyper-prior: Normal or Jeffreys
#
# @return       a sample of p, mu and sigma and latent group membership
#
# @details      Given a sample of plausible values from a mixture of two normal distributions
# this function produces one sample of means and variances of plausible values from their posterior
#
# Straightforward implementation using conjugate priors. See Chapter 6 of 
# Marin, J-M and Robert, Ch, P. (2014). Bayesian essentials with R. 
# 2nd Edition. Springer: New-York. Or Section 6.2.1 of 
# Frühwirth-Schnatter, S. (2006). Finite mixture and Markov switching models. Springer.
#
# Our implementation of Jeffreys prior based on Marsman, M., et al. (2018) An introduction to network 
# psychometrics: Relating Ising network models to item response theory models. 
# Multivariate behavioral research 53.1 (2018): 15-35.
#
update_pv_prior_mixnorm = function (pv, p, mu, sigma, prior.dist=c('normal','Jeffreys')) 
{
  prior.dist = match.arg(prior.dist)
  n = length(pv)
  mean_pv = mean(pv)
  var_pv = var(pv)
  z = rep(0, n)
  nj = c(0,0); 

  if (prior.dist == 'normal')
  {
     ## hyper-prior
    l = rep(1,2)
    v = rep(5,2)

      ## latent group membership
    for (t in 1:n)
    {
      prob = p*dnorm(pv[t], mean = mu, sd = sigma)
      if (sum(prob)==0) prob=c(0.5,0.5)
      z[t] = sample(c(1,2), size = 1, prob=prob)
    }

      ## Means
    for (j in 1:2)
    {
        nj[j]  = sum(z==j) 
        mu[j]  = rnorm(1, mean = (l[j]*mean_pv+nj[j]*mean(pv[z==j]))/(l[j]+nj[j]), 
                           sd = sqrt(sigma[j]^2/(nj[j]+l[j])))
    }
    
      ## Variance
    for (j in 1:2) 
    {
        if (nj[j]>1)
        {
          var_pvj = var(pv[z==j])
          sigma[j] = sqrt(1/rgamma(1, shape = 0.5*(v[j]+nj[j]) ,
                                      rate = 0.5*(var_pv + nj[j]*var_pvj + 
                                               (l[j]*nj[j]/(l[j]+nj[j]))*(mean_pv - mean(pv[z==j]))^2)))
        }
    }
  
      ## membership probabilities
    p = rgamma(2, shape = nj + 1, scale = 1)
    p = p/sum(p)
  }
  
  if (prior.dist == 'Jeffreys')
  {
      ## latent group membership
      prob = p[1] * dnorm(pv, mu[1], sigma[1])
      prob = prob / (prob + (1-p[1]) * dnorm(pv, mu[2], sigma[2]))
      z = rbinom (n, 1, prob) + 1
    
        ## Means
      for (j in 1:2)
      {
        nj[j]  = sum(z==j) 
        mu[j]  = rnorm(1, mean = mean(pv[z==j]), 
                       sd = sigma[j]/sqrt(nj[j]))
      }
        ## Variances
      for (j in 1:2) 
      {
        if (nj[j]>1) sigma[j] = 1/sqrt(rgamma(1, nj[j]/2, sum((pv[z==j]-mu[j])^2)/2 ))
      }
    
        # Membership probabilities
      p[1] = rbeta(1, nj[1] + 0.5, nj[2] + 0.5)
      p[2] = 1-p[1]
  }

  return(list(p = p, mu = mu, sigma = sigma, grp = z))
}

####################### MAIN FUNCTION TO GENERATE PLAUSIBLE VALUES

# Plausible values. 
# @param x                tibble(booklet_id <char or int>, person_id <char>, booklet_score <int>, pop <int>)
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
# tibble(booklet_id <char or int>, person_id <char>, booklet_score <int>, nPV nameless columns with plausible values)
pv = function(x, design, b, a, nPV, from = NULL, by = NULL, prior.dist = c("normal", "mixture"))
{
  prior.dist = match.arg(prior.dist)
  nPop = length(unique(x$pop))
  n_prior_updates = 10
  if (is.null(from)) from = Gibbs.settings$from.pv
  if (is.null(by))   by = Gibbs.settings$step.pv
  
  pb = get_prog_bar(nsteps = n_prior_updates+nPV)
  on.exit({close_prog_bar()})
  

  if (prior.dist == "mixture")
  {
    priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2))
    x$grp = sample(1:2, nrow(x), replace=T, prob=c(0.5,0.5))
  } else
  {
    priors = list(mu=rep(0,nPop), sigma=rep(4,nPop), mu.a=0, sigma.a=1)
  }
  
  if (is.matrix(b))
  {
    which.pv = seq(from,(from-by)*(from>by)+by*nPV,by=by)
    nIter=max(which.pv)
    if (nrow(b)<nIter){
      stop(paste("at least", as.character(nIter), "samples of item parameters needed in function pv"))
    }
    b.step = as.integer(nrow(b)/nIter)
    
    pb$set_nsteps(nIter)
    for(iter in 1:nIter)
    {
      if (prior.dist == "mixture")
      {
        x = x %>% 
          group_by(.data$booklet_id, .data$grp) %>%
          mutate(PVX = pv_recycle(b[iter*b.step,], a, 
                                  design[[.data$booklet_id[1]]]$first, 
                                  design[[.data$booklet_id[1]]]$last, 
                                  .data$booklet_score, 1, 
                                  priors$mu[.data$grp[1]], priors$sigma[.data$grp[1]])) %>%
          ungroup() 
        
        priors = update_pv_prior_mixnorm(x$PVX, priors$p, priors$mu, priors$sigma)
        x$grp = priors$grp
      } else if (prior.dist == "normal") 
      {
        x = x %>% 
          group_by(.data$booklet_id,.data$pop) %>%
          mutate(PVX = pv_recycle(b[iter*b.step,], a, 
                                  design[[.data$booklet_id[1]]]$first, 
                                  design[[.data$booklet_id[1]]]$last, 
                                  .data$booklet_score, 1, 
                                  priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
          ungroup()
        
        if (nPop==1) priors = update_pv_prior(x$PVX, x$pop, priors$mu, priors$sigma)
        if (nPop>1)  priors = update_pv_prior_H(x$PVX,x$pop,priors$mu, priors$sigma, priors$mu.a, priors$sigma.a)
      }
      
      if (iter %in% which.pv)
      {
        colnames(x)[colnames(x)=='PVX'] = paste0('PV', iter)
      }
      pb$tick()
    }
    return(select(x, .data$booklet_id, .data$person_id,.data$booklet_score, matches('PV\\d+')))
    
  }else # if b is not a matrix
  {
    # it is safe to use the ordered pv's for the prior update
    for(iter in 1:n_prior_updates) 
    {
      if (prior.dist == "mixture")
      {
        x = x %>% 
          group_by(.data$booklet_id, .data$grp) %>%
          mutate(PVX = pv_recycle_sorted(b, a, 
                          design[[.data$booklet_id[1]]]$first, 
                          design[[.data$booklet_id[1]]]$last, 
                          .data$booklet_score, 1, 
                          priors$mu[.data$grp[1]], priors$sigma[.data$grp[1]])) %>%
          ungroup() 
        
        priors = update_pv_prior_mixnorm(x$PVX, priors$p, priors$mu, priors$sigma)
        x$grp = priors$grp
      } else if (prior.dist == "normal")
      {
        pv = x %>%
          group_by(.data$booklet_id,.data$pop) %>%
          mutate(pv = pv_recycle_sorted(b, a, 
                      design[[.data$booklet_id[1]]]$first, design[[.data$booklet_id[1]]]$last, 
                      .data$booklet_score, 1L, priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
          ungroup() 

        if (nPop==1) priors = update_pv_prior(pv$pv,pv$pop, priors$mu, priors$sigma)
        if (nPop>1)  priors = update_pv_prior_H(pv$pv,pv$pop,priors$mu, priors$sigma, priors$mu.a, priors$sigma.a)
      }
      pb$tick()
    }
    return(
      x %>% 
          group_by(.data$booklet_id, .data$pop) %>%
          do({
            bkID = .$booklet_id[1]
            popnbr = .data$pop[1]
            out_pv = pv_recycle(b, a, design[[bkID]]$first, design[[bkID]]$last, .$booklet_score, 
                                nPV, priors$mu[popnbr], priors$sigma[popnbr])
            data.frame(.$person_id, .$booklet_score, as.data.frame(out_pv), stringsAsFactors = FALSE)
          }) %>%
          ungroup() %>%
          select(-.data$pop)
    )
  }
}


# borrowed the following from mvtnorm source for the time being
# so we don't have to import a whole package
# will make a cpp implementation for next version so this can be removed
rmvnorm = function(n,mean,sigma)
{
  ev = eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  retval = matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% R
  retval = sweep(retval, 2, mean, "+")
  colnames(retval) = names(mean)
  retval
}


update_MVNprior = function(pvs,Sigma)
{
  m_pv = colMeans(pvs)
  nP = nrow(pvs)
  mu = rmvnorm(1, mean=m_pv,sigma=Sigma/nP)
  
  S=(t(pvs)-m_pv)%*%t(t(pvs)-m_pv)
  Sigma = solve(rWishart(1,nP-1,solve(S))[,,1]) 
  return(list(mu = mu, Sigma = Sigma))
}
