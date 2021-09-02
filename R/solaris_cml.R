


swap_sol = function(x) c(x[2],x[1])


elsym_sol = function(b,a,first,last)
{
  n=length(first)
  Msc=0L
  gg=matrix(0,sum(a[last])+1,2)
  gg[1,1]=1
  col=1:2 
  for (i in 1:n)
  {
    #for (s in 0:Msc) gg[s+1,col[2]]=0
    gg[,col[2]]=0
    for (s in 0:Msc)
    {
      for (j in first[i]:last[i])
      {
        gg[s+1+a[j],col[2]]=gg[s+1+a[j],col[2]]+gg[s+1,col[1]]*b[j]
      }
    }
    Msc = Msc+a[last[i]];
    col = swap_sol(col)
  }

  return(gg[,col[1]])
}

E.STEP_sol = function(b, a , first, last, scoretab)
{
  n = length(first)
  ms= length(scoretab)-1
  expect = rep(0,length(b)) 
  g = elsym_sol(b,a,first,last)
  for (i in 1:n)
  {
    gi = elsym_sol(b,a,first[-i],last[-i])
    for (j in first[i]:last[i])
    {
      for (s in 1L + a[j]:ms)
      {
        #s = s+1
        if ((g[s]>0)&&((s-a[j])<ms)) 
          expect[j] = expect[j] + scoretab[s]*(gi[s-a[j]]*b[j]/g[s])
      }
    }
  }
  return(expect)
}

H.STEP_sol = function(b, a, first, last, scoretab)
{
  nI=length(last)
  ms=length(scoretab)-1
  H=matrix(0,length(a),length(a))
  
  g=elsym_sol(b,a,first,last)
  for (item in 1:nI)
  {
    gi=elsym_sol(b,a,first[-item],last[-item])
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
          gk=elsym_sol(b,a,first[-k],last[-k])
          gik=elsym_sol(b,a,first[-c(item,k)],last[-c(item,k)])
          for (l in (first[k]+1):last[k])
          {
            for (s in 1:ms)
            {
              if (g[s]>0)
              {
                if ((s>(a[j]+a[l]))&&((s-a[j]-a[l])<=length(gik))){
                  H[l,j] = H[l,j] + scoretab[s]*(gik[s-a[j]-a[l]])*((b[j]*b[l])/g[s]) 
                }
                if ((s>a[j])&&(s>a[l])) {
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

est_lambda_sol <- function(b, a, first, last, scoretab)
{
  ifelse(scoretab>0, scoretab/(elsym_sol(b,a,first,last)*sum(scoretab)), NA) 
}  

#calibrate_Bayes(scoretab, design, sufI, a, first, last,  nIter, fixed_b=fixed_b)

# to do: aanroep aanpassen, list of booklets wordt niet meer gebruikt
# calibrate_CML_sol = function(booklet, sufI, a, first, last, nIter, fixed_b=NULL) {
calibrate_CML_sol = function(scoretab, design, sufI, a, first, last, nIter, fixed_b=NULL)
{
  #nb = n_distinct(scoretab$booklet_id)
  ni = length(first)
  max_nr_iter = 30
  
  EsufI = rep(0,length(sufI))
  
  bk_design = split(design, design$booklet_id,drop=TRUE)
  bk_scoretab = split(scoretab, scoretab$booklet_id,drop=TRUE)
  
  pb = get_prog_bar()
  on.exit({pb$close()})
  
  if (is.null(fixed_b)) # if no fixed parameters
  {
    nn= sum(sufI)
    b = rep(1,length(a))
    ## Implicit Equations  ###
    converged=FALSE
    iter=0

    while ((!converged)&&(iter<=nIter))
    {
      iter=iter+1
      EsufI[] = 0
      for(bk in bk_design)
      {
        EsufI = EsufI + E.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
      }
      b = b*sufI/EsufI
      converged=(max(abs(sufI-EsufI))/nn<1e-04)
      if(is.na(converged))
      {
        return(calibrate_Bayes(scoretab, design, sufI, a, first, last,  nIter, fixed_b=fixed_b))
      }
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
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1
    while ((!converged)&&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI[] = 0
      H[] = 0
      for(bk in bk_design)
      {
        EsufI = EsufI + E.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
        H     = H     + H.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
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

      b = try(b*exp(solve(H*scale,sufI-EsufI)))
      converged=(max(abs(EsufI-sufI))/nn<1e-10)
      
      if(inherits(b,'try-error') || is.na(converged))
      {
        return(calibrate_Bayes(scoretab, design, sufI, a, first, last,  nIter, fixed_b=fixed_b))
      }
      
      pb$tick()
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

    m=sum(scoretab$N)
    
    nn=m*ni_free
    
    converged=FALSE
    iter=0
    pb = txtProgressBar(min=0, max=nIter)
    while ((!converged)&&(iter<=nIter))
    {
      iter=iter+1
      EsufI[] = 0
      for(bk in bk_design)
      {
        EsufI = EsufI + E.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
      }
      b[update_set] = b[update_set]*sufI[update_set]/EsufI[update_set]
      converged=(max(abs(sufI[update_set]-EsufI[update_set]))/nn<1e-04)
      if(is.na(converged))
      {
        return(calibrate_Bayes(scoretab, design, sufI, a, first, last,  nIter, fixed_b=fixed_b))
      }
      pb$tick()
    }
    ie_iter=iter
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1
    while ((!converged)&&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI[] = 0
      H[] = 0
      for(bk in bk_design)
      {
        EsufI = EsufI + E.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
        H     = H     + H.STEP_sol(b,a,bk$first,bk$last,bk_scoretab[[bk$booklet_id[1]]]$N)
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
      b = try(b*exp(solve(H*scale,sufI-EsufI)))
      converged=(max(abs(EsufI[update_set]-sufI[update_set]))/nn<1e-10)
      
      if(inherits(b,'try-error') || is.na(converged))
      {
        return(calibrate_Bayes(scoretab, design, sufI, a, first, last,  nIter, fixed_b=fixed_b))
      }
      
      pb$tick()
      scale=1
    }
    close(pb)
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
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
    bk_design, bk_scoretab,
    SIMPLIFY=FALSE, USE.NAMES=FALSE) %>%
    bind_rows()
  
  return(list(b=b, H=H, beta=report$beta, acov.beta=report$cov.beta,
              lambda=lx, n_iter=iter, nr_iter = nr_iter, ie_iter=ie_iter))
}


