
# returns a simplified parms object mainly useful for person ability estimation
# parms can be data frame or a parms object
#
# returns list: a, b, design (including first and last), items(including first and last), method (CML or Bayes)
#
# draw can be a number
# the default 'sample' returns the complete matrix, i.e. the sample is the census in the original order
# any subsequent 'sampling' must be done in the caller function itself

simplify_parms = function(parms, design=NULL, draw = c('sample','average'), by_chain=FALSE)
{
  check_df(design, 'item_id', nullable=TRUE)
  chain_index = NULL
  if(!is.numeric(draw))
  {
    draw = match.arg(draw)
  } else
  {
    draw = as.integer(draw)
    check_num(draw,'integer',.length=1)
  }
  
  if(!is.null(design))
  {
    design = ungroup(design)
    design$item_id = as.character(design$item_id)
    if(!('booklet_id' %in% colnames(design)))
      design$booklet_id='all_items'
    
    design = design[,c('booklet_id','item_id')]
    if(n_distinct(design$item_id,design$booklet_id) < nrow(design))
    {
      warning('Duplicated items are removed in one or more booklets in design')
      design = distinct(design)
    }
  }
  
  if(inherits(parms, 'mst_enorm') && !'design' %in% names(parms$inputs))
  {
    if(is.null(design))
      design = lapply(parms$inputs$bkList,function(bk) tibble(item_id=as.character(bk$items))) |>
        bind_rows(.id='booklet_id')
    
    parms = coef(parms)
  }
  
  if(inherits(parms,'enorm') || inherits(parms,'prms'))
  {
    method = parms$inputs$method
    a = parms$inputs$ssIS$item_score
    if(method=="CML"){
      b = parms$est$b
    } else if(is.numeric(draw)) 
    {
      if(draw>nrow(parms$est$b)) 
        stop_(sprintf('Draw %i is larger than the number of item parameter draws (%i)', draw,nrow(parms$est$b)))
      b = parms$est$b[draw,]	
    } else if(draw == 'average')
    {
      if(by_chain)
      {
        s = parms$est$chain_start
        chain_index = mapply(s,c(s[-1]-1,nrow(b)), FUN=`:`,SIMPLIFY=FALSE)
        b = do.call(rbind, lapply(chain_index, function(indx) colMeans(b[indx,])))
      } else
      {
        b = colMeans(parms$est$b)  
      }
    } else 
    {
      b = parms$est$b
    } 
    
    fl = parms$inputs$ssI[,c('item_id','first','last')]
    
    if(is.null(design))
    {
      design = parms$inputs$design
    } else
    {
      old_itm = design$item_id
      design$item_id = factor(as.character(design$item_id), levels=parms$inputs$ssI$item_id)
      if(anyNA(design$item_id))
      {
        message("the following item_id's in the data or design are not present in your parameters")
        print(setdiff(as.character(old_itm), levels(parms$inputs$ssI$item_id)))
        stop('items without parameters') 
      }
      design = inner_join(design, fl, by='item_id')
    }
  } else
  {
    method='CML'
    parms$item_id = as.character(parms$item_id)
    # to do: Assumption that caller is not interested in other items when a design is specified
    # but this is not done when parms is not a df
    # all such uses atmo are indeed correct in both cases and it does speed things up... maybe use an extra argument?
    if(!is.null(design)) parms = semi_join(parms,design,by='item_id')
    
    parms = transform.df.parms(parms,'b')
    parms$item_id = factor(parms$item_id)
    # transform.df result is guaranteed sorted
    a = parms$item_score
    b = parms$b
    
    if(n_distinct(parms$item_id)==1)
    {
      fl = tibble(item_id=parms$item_id[1],first=1L,last=nrow(parms))
    } else
    {
      fl = parms |>
        mutate(rn=row_number()) |>
        group_by(.data$item_id) |>
        summarise(first=as.integer(min(.data$rn)), last=as.integer(max(.data$rn))) |>
        ungroup() |>
        arrange(.data$first)
    }
    
    if(is.null(design))
    {
      design = mutate(fl, booklet_id='all_items')
    } else
    {
      old_itm = design$item_id
      design$item_id = ffactor(as.character(design$item_id), levels=levels(parms$item_id))
      if(anyNA(design$item_id))
      {
        message("the following item_id's in the data or design are not present in your parameters")
        print(setdiff(as.character(old_itm), levels(parms$item_id)))
        stop('items without parameters') 
      }
      design = inner_join(design, fl, by='item_id')
    }
  }
  
  fl$first0 = fl$first - 1L
  fl$last0 = fl$last - 1L
  design$first0 = design$first - 1L
  design$last0 = design$last - 1L
  
  design = arrange(design, .data$booklet_id,.data$first)
  
  booklets = design |>
    group_by(.data$booklet_id) |>
    summarise(nit = n(), max_score = sum(a[.data$last]) ) |>
    ungroup() |>
    arrange(.data$booklet_id)

  list(a=a, b=b, booklets=booklets, design=design, items = fl, method=method, chain_index=chain_index)
}

# item_id is returned as a factor with the same levels as deign$item_id
fixed_param2df = function(fixed_params, design)
{
  if(!is.null(fixed_params))
  {
    if(inherits(fixed_params,'enorm') || inherits(fixed_params,'prms'))
    {
      if(inherits(fixed_params,"mst_enorm"))
      {
        fixed_params$inputs = fixed_params$mst_inputs
        fixed_params$est$b = fixed_params$mst_est$b
        fixed_params$est$a = fixed_params$mst_est$a
      }
      
      if (fixed_params$inputs$method!="CML")
        message("Posterior means are taken as values for fixed parameters")
      
      fixed_params = fixed_params$inputs$ssIS
      fixed_params$b = if.else(fixed_params$inputs$method=="CML", fixed_params$est$b, colMeans(fixed_params$est$b))
      
    } else
    {
      # transform the fixed params to the b parametrization dexter uses internally
      fixed_params = transform.df.parms(fixed_params, out.format = 'b') 
    }  
    
    # avoid factor warnings and reduce
    fixed_params$item_id = factor(as.character(fixed_params$item_id), levels=levels(design$item_id))
    
    fixed_params = filter(fixed_params,!is.na(.data$item_id))  |>
      select('item_id','item_score','b') |>
      arrange(.data$item_id,.data$item_score)

  }
  fixed_params
}



transform.df.parms = function(parms.df, out.format = c('b','beta','eta'))
{
  out.format = match.arg(out.format)
  colnames(parms.df) = tolower(colnames(parms.df))
  parms.df = ungroup(parms.df)
  
  if('delta' %in% colnames(parms.df))
    parms.df = rename(parms.df, beta = 'delta')
  in.format = intersect(colnames(parms.df), c('b','beta','eta'))
  
  if(length(in.format) == 0)
    stop('parameters must contain  one of: b, beta, eta')
  
  if(length(in.format)>1)
  {
    in.format = in.format[1]
    message(paste0("Using '",in.format,"' as input parameter"))
  }
  
  if(!all(c('item_id','item_score') %in% colnames(parms.df)))
    stop('parameters must contain the columns: item_id, item_score')
  
  if(any(parms.df$item_score%%1 > 0))
    stop("column 'item_score' must be integer valued")
  
  if(n_distinct(parms.df$item_id, parms.df$item_score) < nrow(parms.df))
    stop('multiple parameters supplied for the same item and score')
  
  
  if(any(parms.df$item_score <= 0))
    stop("Items scores of 0 or less are not supported")
  
  if(in.format == 'b' && any(parms.df$b <= 0))
    stop("A 'b' parameter cannot be negative, perhaps you meant to include a 'beta' parameter?")
  
  parms.df$item_id = as.character(parms.df$item_id)
  parms.df$item_score = as.integer(parms.df$item_score)
  
  if(in.format!=out.format)
  {
    if(in.format == 'beta' && out.format == 'b')
    {
      if(n_distinct(parms.df$item_id) == 1) 
      {
        parms.df$b = item_beta2b(parms.df$item_score,parms.df$beta)
        return(parms.df)
      }
      parms.df = parms.df |>
        group_by(.data$item_id) |>
        arrange(.data$item_score) |>
        mutate(b = item_beta2b(.data$item_score,.data$beta)) |>
        ungroup()
    } else if(in.format == 'beta' && out.format == 'eta')
    {
      parms.df = parms.df |>
        group_by(.data$item_id) |>
        arrange(.data$item_score) |>
        mutate(b = item_beta2eta(.data$item_score,.data$beta)) |>
        ungroup()
    } else if(in.format== 'eta' && out.format == 'b')
    {
      parms.df$b = exp(-parms.df$eta)
    } else
    {
      stop_(sprintf('transformation %s -> %s not implemented', in.format, out.format))
    }
  } 
  # must be sorted!
  arrange(parms.df,.data$item_id, .data$item_score)
}





# parametrisations --------------------------------------------------------

# a and beta must be sorted by a

item_beta2b = function(a,beta)
{
  exp(-cumsum(beta * (a - c(0L,a[-length(a)]))))
}


item_beta2eta = function(a,beta)
{
  cumsum(beta * (a - c(0L,a[-length(a)])))
}


beta2b = function(a,beta,first,last)
{
  indx = mapply(first,last,FUN=':',SIMPLIFY=FALSE)
  b = double(length(beta))
  for(i in indx)
  {
    b[i] = item_beta2b(a[i],beta[i])
  }
  b
}


# from dexter
makeD = function(a,first,last)
{
  k = length(a)
  D = matrix(0,k,k)
  tel=1
  for (i in 1:length(first))
  {
    for (j in 1:(last[i]-first[i]+1))
    {
      if (j==1){
        D[tel,tel]=-1/a[tel]
      }else
      {
        D[tel,tel-1]=-1/(a[tel-1]-a[tel])
        D[tel,tel]=1/(a[tel-1]-a[tel])        
      }
      tel=tel+1
    }
  }
  return(D)
}

# copied from previous dexter version, just removed the part that removes the 0 cat
toOPLM = function(a, b, first, last, H=NULL, fixed_b=NULL, lambda=NULL, method=c('CML','Bayes'))
{
  b_rn = b
  a_org = a

  logb = log(b)
  cov.beta = NULL
  method=match.arg(method)

  ### Bayesian: b is a matrix
  if(method=='Bayes')
  {
    k=ncol(b)
    beta=b
    for (r in 1:nrow(b))
    {
      for (i in 1:length(first))
      {
        beta[r,first[i]]=-logb[r,first[i]]/a[first[i]]
        if ((last[i]-first[i])>0)
        {
          tmp=(logb[r,(first[i]+1):last[i]]-logb[r,first[i]:(last[i]-1)])
          tmp=tmp/(a[first[i]:(last[i]-1)]-a[(first[i]+1):last[i]])
          beta[r,(first[i]+1):last[i]]=tmp
        }
      }
      if (is.null(fixed_b)){
        mean_beta = mean(beta[r,])
        b_rn[r,] = b_rn[r,]*exp(mean_beta*a_org)
        beta[r,] = beta[r,] - mean_beta ## mean center
      }
    }
    if (nrow(beta)>2) cov.beta=cov(beta)
  }else if(method=='CML')
  {                                       ### CML; b is a vector
    DD = makeD(a,first,last)
    if (is.null(fixed_b))
    {
      beta = DD%*%logb
      b_rn = b_rn*exp(mean(beta)*a_org) # re-normalize b such that it corresponds to beta
      k  = length(b)
      CC = matrix(-1/k,k,k)
      diag(CC)=(k-1)/k
      beta = CC%*%beta
      if (!is.null(H))
      {
        A  = CC%*%DD
        cov.beta = solve(H[-1,-1])
        cov.beta = A[,-1]%*%cov.beta%*%t(A[,-1])
      }
    }else # if there are fixed parameters we do not (re-)normalize
    {
      beta = DD%*%logb
      if (!is.null(H))
      {
        cov.beta  = solve(H[!fixed_b,!fixed_b])
        cov.beta  = DD[,!fixed_b]%*%cov.beta%*%t(DD[,!fixed_b]) #this brings the fixed params back
      }
    }
  }  
  return(list(beta=beta, cov.beta=cov.beta, a=a, b_renorm = b_rn, first=first, last=last))
}


