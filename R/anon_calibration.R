

elsym = function(b,a,first,last,i1=0L,i2=0L)
{
  n = length(first)
  ms = as.integer(sum(a[last]))
  g = double(ms+1L)
  
  ElSym_C(b, a, as.integer(first-1L), as.integer(last-1L), 
          as.integer(i1-1L), as.integer(i2-1L), n, ms, g)

  return(g)
}


# returns log likelihood, or vector of log likelihoods if bayes
logL = function(parms, mean_gibbs=FALSE)
{
  a = parms$inputs$ssIS$item_score
  b = parms$est$b
  if(is.matrix(b) && mean_gibbs)
    b=colMeans(b)
  
  scoretab = split(parms$inputs$scoretab,parms$inputs$scoretab$booklet_id,drop=FALSE)
  design = split(parms$inputs$design, parms$inputs$design$booklet_id, drop=TRUE)
  
  llb = function(b)
  {
    ll_rm = sum(parms$inputs$ssIS$sufI * log(b))
    
    lgRM = sapply(scoretab,
           function(stb)
           {
             d = design[[stb$booklet_id[1]]]
             log(elsym(b, a, d$first, d$last)) * stb$N
           }) %>%
      unlist() %>%
      sum()
  
    ll_rm-lgRM
  }
  if(is.matrix(b))
  {
    apply(b,1, llb)
  } else
  {
    llb(b)
  }
}


######################################################################
# Estimation of Rasch and Interaction Model 
#######################################################################
# Using Newton-Raphson per item
# Currently with elsym and scale
# @returns:
#      bRM:    Parameter estimates of the Rasch model
#      bIM:    Parameter estimates of the Interaction model
#      cIM:    Estimate of (exp-)interaction parameter
#      cRM:    Interaction parameters under the Rasch model: all equal to 1
#      HRM:    Block diag. Asymptotic var-covar matrix of item parameters under RM
#     se.sigma:   Standard error of interaction parameter
# fit.stat: log(cIM)/se.sigma. Wald statistic normally distributed under Rasch model
#########################################################################
# future~to~d0: switch to ittotmat mean if overflow, use full Hessian

EstIM = function(first,last, nCat, a, sufI, sufC, scoretab, regs=FALSE) {
  
  C = rep(1:length(first), nCat)
  
  m=sum(scoretab) ##
  nI=length(last)
  b=rep(0,length(sufI))
  ic=rep(1,nI)
  var.ic = double(nI)
  #HRM=matrix(0,length(b),length(b))
  HIM = vector("list", nI)
  
  # Identification
  b[sufI>0]=1
  upd_set=vector("list",nI)
  for (i in 1:nI)
  {
    upd_set[[i]]=which(sufI[first[i]:last[i]]>0)
    upd_set[[i]]=upd_set[[i]][-1]
    upd_set[[i]]=(first[i]:last[i])[upd_set[[i]]]
  }
  
  ps = possible_scores(a, first, last)
  # see actual scores, do: (1:length(ps)-1L)[as.logical(ps)]
 
  converged=2
  scale=2

  while(converged>0.01)
  {
    converged=-1
    pi_mat = ittotmat0(b,ic[C],a,first,last, ps)
    
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
        converged=pmax(converged,max(abs(E))/m) #
        #HRM[upd_set[[i]],upd_set[[i]]]=H
      }
    }
    if (converged<1) scale=1
  }
  
  bRM=b
  cRM=ic
  
  ## IM
  converged=2
  scale=2
  first_iter = TRUE
  while(converged>0.001)
  {
    converged=-1
    pi_mat = ittotmat0(b,ic[C],a,first,last, ps) 

    if(first_iter)
    {
      first_iter = FALSE
      ctrRM = pi_mat
    }
    for (i in 1:nI)
    {
      # gradient and hessian for thresholds of item i
      if (length(upd_set[[i]])>0)
      {
        pi = pi_mat[upd_set[[i]], , drop=FALSE]
        E = sufI[upd_set[[i]]] - pi %*% scoretab
        H = -pi %*% tcrossprod(diag(scoretab), pi)
        diag(H) = pi %*% scoretab + diag(H)
        
        # gradient and hessian for interaction parameter
        ncol_pi = ncol(pi)
        nrow_pi = nrow(pi)
        s_range = 0:(ncol_pi-1)
        
        E = c(E, sufC[i])
        H = cbind(H, rep.int(0, nrow(H)))
        H = rbind(H, rep.int(0, ncol(H)))
        k = 1
        e0 = 0; e1 = 0
        f = matrix(0, nrow_pi, ncol_pi)
        h = 0
        for (j in upd_set[[i]])
        {
          E[length(E)]=E[length(E)]-a[j]*sum(s_range*scoretab*pi[k,])
          e0=e0+a[j]*pi[k,]
          e1=e1+a[j]^2*pi[k,]
          f[k,]=a[j]*s_range*pi[k,]
          h=h+f[k,]
          k=k+1
        }
        H[nrow(H),nrow(H)]=sum(s_range^2*(e1-e0^2)*scoretab)
        for (k in 1:nrow(f))
        {
          H[k,nrow(H)]=sum((f[k,]-pi[k,]*h)*scoretab)
          H[nrow(H),k]=H[k,nrow(H)]
        }
        
        
        # NR update for parameters of item i
        update = solve(H*scale,E)
        b[upd_set[[i]]] = b[upd_set[[i]]]*exp(update[-length(update)])
        ic[i] = ic[i]*exp(update[length(update)])
        HIM[[i]] = H
        var.ic[i] = solve(H)[nrow(H),nrow(H)]
        converged = pmax(converged,max(abs(E))/m)
      }
    }
    if (converged<1) scale=1
  }
  
  sigma = log(ic)
  sigma = sigma - mean(sigma)
  ic = exp(sigma)
  se.sigma = sqrt(var.ic)
  fit.stats = sigma/se.sigma
  
  out = list(bRM=bRM,cRM=cRM,bIM=b,cIM=ic,se.sigma=se.sigma,HIM=HIM, fit.stats=fit.stats, possible_scores = (1:length(ps)-1L)[as.logical(ps)])
  if(regs)
  {
    out$ctrRM = ctrRM
    out$ctrIM = ittotmat0(b,ic[C],a,first,last, ps) 
  }
  # ###### test ####### #
  #HH=H_IM(a, b, ic, first, last, sufI, sufC, scoretab, method="full")
  #HR=solve(HH$Hessian)
  #plot(fit.stats, log(ic)/sqrt(diag(HR)[last])); abline(0,1,lty=2)
    # H_full = matrix(0, length(b), length(b))
    # Grad_full = double(length(b))
    # pi_s_full = matrix(0, length(b), length(scoretab))
    # H_im(a,b,ic[C], as.integer(first-1L), as.integer(last-1L), sufI, sufC, scoretab, H_full, Grad_full, pi_s_full, FALSE)
  # # # 
    # plot(diag(solve(H_full))[last], diag(HR)[last]); abline(0,1,lty=2)
    # H_diag = matrix(0, length(b), length(b))
    # Grad_diag = double(length(b))
    # pi_s_diag = matrix(0, length(b), length(scoretab))
    # H_im(a,b,ic[C], as.integer(first-1L), as.integer(last-1L), sufI, sufC, scoretab, H_diag, Grad_diag, pi_s_diag, TRUE)
    # plot(diag(solve(H_diag))[last], diag(HR)[last]); abline(0,1,lty=2)
  #  for (i in 1:24){print(max(abs(H_diag[first[i]:last[i],first[i]:last[i]]-HIM[[i]])))}
  # browser()
  # ###### end test ####### #
  out
}






# a, first, last, sufI: all ordered by item_id, item_score. include 0 score
# scoretab: data.frame(booklet_id, booklet_score, N), including scores with N=0, ordered by booklet_id, booklet_score
# design: data.frame(booklet_id, item_id, first, last), ordered by booklet_id, first
# fixed_b:  NUll: nothing fixed; 
#           vector length of b with NA values for parameters that need to be estimated, b values for fixed 

E_bkl = function(..., use_mean = FALSE)
{
  if(use_mean)
  {
    E_booklets_mean(...)[,1,drop=TRUE]
  } else
  {
    E_booklets(...)[,1,drop=TRUE]
  }
}

NR_bkl = function(..., use_mean = FALSE)
{
  if(use_mean)
  {
    NR_booklets_mean(...)
  } else
  {
    NR_booklets(...)
  }
  invisible(NULL)
}


calibrate_CML = function(scoretab, design, sufI, a, first, last, nIter=1000, fixed_b=NULL) 
{
  pb = get_prog_bar()
  on.exit({pb$close()})

  use_mean = FALSE

  # bookkeeping, make counts for C functions
  # items per booklet
  nib = design %>%
    count(.data$booklet_id) %>%
    arrange(.data$booklet_id) %>%
    pull(.data$n)
  
  # number of score per booklet (max_score + 1)
  nscore = scoretab %>%
    count(.data$booklet_id) %>%
    arrange(.data$booklet_id) %>%
    pull(.data$n)
  
  # THESE ARE C INDEXED!
  # one vector for all
  bfirst = as.integer(design$first - 1L)
  blast = as.integer(design$last - 1L)
  
  
  nb = length(nib)
  ni = length(first)
  max_nr_iter = 30
  
  max_par_bk = design %>% 
    group_by(.data$booklet_id) %>%
    summarise(npar=sum(.data$last - .data$first+1L)) %>%
    pull(.data$npar) %>%
    max() %>%
    as.integer()
  
  # end bookkeeping

  if (is.null(fixed_b)) # if no fixed parameters
  {
    nn= sum(sufI)
    b = rep(1,length(a))
    ## Implicit Equations  ###
    converged=FALSE
    iter=0
    
    while (!converged && iter<=nIter)
    {
      iter=iter+1
      EsufI = E_bkl(b, a, bfirst, blast, scoretab$N, nscore, nib, use_mean=use_mean)

      converged=(max(abs(sufI-EsufI))/nn < 1e-04)
      if(is.na(converged))
      {
        if(use_mean)
          stop('problem cannot be computed')
        # continue with the next iteration using elsym_mean
        use_mean = TRUE
        converged = FALSE
        next
      } 
      b = b*sufI/EsufI
      pb$tick()
    }
    
    ie_iter=iter
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
    ### identification ###
    # within items
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    # between items
    ref_cat=2
    b[-first] = b[-first]/(b[ref_cat]^(a[-first]/a[ref_cat]))

    ###  NR  ###
    H = matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1 #to~do: ask timo, scale does not do anything at all
    while (!converged && nr_iter<max_nr_iter)
    {
      iter=iter+1
      nr_iter=nr_iter+1

      # sets EsufI and H to 0 and updates in place 
      NR_bkl(b, a, bfirst, blast, scoretab$N, nscore, nib, max_par_bk, EsufI, H, use_mean=use_mean)
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
      
      converged=(max(abs(EsufI-sufI))/nn<1e-10)
      nb = try(b*exp(solve(H*scale,sufI-EsufI)))
      
      if(inherits(nb,'try-error') || is.na(converged))
      {
        if(use_mean)
          stop('problem cannot be computed')
        # continue with the next iteration using elsym_mean
        use_mean = TRUE
        converged = FALSE
        next 
      }
      b = nb
      pb$tick()
      if (nr_iter==2) scale=1
    }
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  }else  ### if fixed parameters
  {
    fixed_set=which(!is.na(fixed_b))
    update_set=which(is.na(fixed_b))
    b=fixed_b
    ni_free=sum(is.na(fixed_b[last]))
    b[update_set]=1
    nn = ni_free * sum(scoretab$N)
    
    converged=FALSE
    iter=0
    while ( !converged && iter<=nIter)
    {
      iter=iter+1
      EsufI = E_bkl(b, a, bfirst, blast, scoretab$N, nscore, nib, use_mean=use_mean)

      converged=(max(abs(sufI[update_set]-EsufI[update_set]))/nn<1e-04)
      
      if(is.na(converged))
      {
        if(use_mean)
          stop('problem cannot be computed')
        # continue with the next iteration using elsym_mean
        use_mean = TRUE
        converged = FALSE
        next
      } 
      
      b[update_set] = b[update_set]*sufI[update_set]/EsufI[update_set]
      pb$tick()
    }
    
    ie_iter=iter
    if (!converged) warning(paste('Implicit Equations not Converged in',nIter,"iterations"))
    
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1
    while (!converged && nr_iter<max_nr_iter)
    {
      iter=iter+1
      nr_iter=nr_iter+1
      
      # sets EsufI and H to 0 and updates in place 
      NR_bkl(b, a, bfirst, blast, scoretab$N, nscore, nib, max_par_bk, EsufI, H, use_mean=use_mean)
      
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
      
      converged=(max(abs(EsufI[update_set]-sufI[update_set]))/nn<1e-10)
      nb = try(b*exp(solve(H*scale,sufI-EsufI)))
      
      if(inherits(nb,'try-error') || is.na(converged))
      {
        if(use_mean)
          stop('problem cannot be computed')
        # continue with the next iteration using elsym_mean
        use_mean = TRUE
        converged = FALSE
        next 
      }
      b = nb
      
      pb$tick()
      scale=1
    }

    if (!converged) warning(paste('Newton-Raphson not Converged in',nr_iter,"iterations"))
  }

  report = toOPLM(a, b, first, last, H=H, fixed_b=fixed_b)
  b = report$b_renorm
  
  lx = mapply(
      function(bdes, bsct)
      {
        tibble(booklet_id = bdes$booklet_id[1],
               booklet_score = bsct$booklet_score,
               lambda = ifelse(bsct$N>0, bsct$N/(elsym(b,a,bdes$first,bdes$last)), NA_real_)) #*sum(bsct$N)
      },
      split(design, design$booklet_id),
      split(scoretab, scoretab$booklet_id), 
      SIMPLIFY=FALSE, USE.NAMES=FALSE) %>%
    bind_rows()
  
  return(list(b=b, H=H, beta=report$beta, acov.beta=report$cov.beta, 
              lambda=lx, n_iter=iter, nr_iter = nr_iter, ie_iter=ie_iter))
}



# scoretab: data.frame: booklet_id, booklet_score, N; ordered by booklet_id, booklet_score; N=0 included
# design: data.frame: booklet_id, item_id, first, last; ordered by booklet_id, first
# fixed_b: NULL or vector of length(b) with NA's for free parameters
# TO DO: At this moment b and lambda are not consistent. b and delta are. We must recalculate the 
# lambda' s using the renormalized b to solve this.
calibrate_Bayes = function(scoretab, design, sufI, a, first, last,  nIter, fixed_b=NULL, 
                           from = 25L, step = 5L,
                           start_b=NULL)
{
  # decide on starting values
  if(!is.null(start_b))
  {
    b = start_b
  } else
  {
    b = calibrate_CML(scoretab, design, sufI, a, first, last, fixed_b = fixed_b)$b 
  }
  prior_eta = 0.5
  prior_rho = 0.5
  
  pb = get_prog_bar(nsteps = nIter)
  on.exit({pb$close()})
  
  # bookkeeping: make counts and indexes for C function
  design = design %>%
    mutate(bn = dense_rank(.data$booklet_id) - 1L) %>%
    group_by(.data$bn) %>%
    mutate(inr = dense_rank(.data$first) - 1L) %>%
    ungroup() %>%
    arrange(.data$first, .data$bn)
  
  bi = design$inr
  ib = design$bn
  
  nbi = design %>% 
    count(.data$first) %>% 
    arrange(.data$first) %>% 
    pull(.data$n)
  
  nib = design %>% 
    count(.data$bn) %>% 
    arrange(.data$bn) %>% 
    pull(.data$n)
  
  design = arrange(design,.data$bn)
  
  bfirst = as.integer(design$first -1L)
  blast = as.integer(design$last -1L)
  
  bmax = design %>%
    group_by(.data$bn) %>%
    summarise(max_score = sum(a[.data$last])) %>%
    ungroup() %>%
    pull(.data$max_score)

  m = scoretab %>%
    group_by(.data$booklet_id) %>%
    summarize(m=sum(.data$N)) %>%
    ungroup() %>%
    arrange(.data$booklet_id) %>%
    pull(.data$m)
  
  
  fixed_b_vec = fixed_b
  if(is.null(fixed_b))
    fixed_b_vec = rep(NA_real_, length(b))
  
  
  nchains = ncores = get_ncores(desired = min(ceiling(nIter/4),128L), maintain_free = 1L)
  b = matrix(rep(b,nchains), ncol=nchains)
  
  #end bookkeeping
  
  out = calibrate_Bayes_chains(as.integer(a), as.integer(first-1L), as.integer(last-1L),
                            ib, bi, nbi, nib, bfirst, blast, bmax, m,
                            sufI, scoretab$N, b, fixed_b_vec, 
                            from, step, 
                            as.integer(nIter), pb$cpp_prog_init(), ncores,
                            prior_eta, prior_rho)


  report = toOPLM(a, out$b, first, last, H=NULL,fixed_b=fixed_b)
  
  colnames(out$lambda) = paste(scoretab$booklet_id, scoretab$booklet_score, sep='-')
  
  return(list(a=a, b=report$b_renorm, lambda=out$lambda, beta=report$beta)) #use b=out$b for consistency with lambda
}

