
# returns a simplified parms object mainly useful for person ability estimation
# parms can be data frame or a parms object
#
# returns list: a, b, design (including first and last), items(including first and last), method (CML or Bayes)
#
# draw can be a number
# the default 'sample' returns the complete matrix, i.e. the sample is the census in the original order
# any subsequent 'sampling' must be done in the caller function itself
# zero indexed refers to c-style indexes. (the 0 score parameter is always included either way)

simplify_parms = function(parms, design=NULL, 
                           draw = c('sample','average'),
                           zero_indexed=FALSE)
{
  check_df(design, 'item_id', nullable=TRUE)
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
    if(!('booklet_id' %in% colnames(design)))
      design$booklet_id='all_items'
    design = design[,c('booklet_id','item_id')]
  }
  
  if(inherits(parms, 'mst_enorm') && !'design' %in% names(parms$inputs))
  {
    if(is.null(design))
      design = lapply(parms$inputs$bkList,function(bk) tibble(item_id=bk$items)) %>%
        bind_rows(.id='booklet_id')
    
    parms = coef(parms)
  }
  
  if(inherits(parms,'prms'))
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
      b = colMeans(parms$est$b)  
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
      design$item_id = ffactor(as.character(design$item_id), levels=parms$inputs$ssI$item_id)
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
    parms = transform.df.parms(parms,'b',include.zero=TRUE)
    parms$item_id = ffactor(as.character(parms$item_id))
    a = parms$item_score
    b = parms$b
    
    fl = parms %>%
      arrange(.data$item_id) %>%
      mutate(rn=row_number()) %>%
      group_by(.data$item_id) %>%
      summarise(first=as.integer(min(.data$rn)), last=as.integer(max(.data$rn))) %>%
      ungroup() 
    
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
  if(zero_indexed)
  {
    fl$first = fl$first - 1L
    fl$last = fl$last - 1L
    design$first = design$first - 1L
    design$last = design$last - 1L
  }
  
  list(a=a, b=b, design=design, items = fl, method=method)
}




# parms.df 
# data.frame with columns item_id, item_score, and one of b, eta, beta/delta
# it is assumed that parms.df comes from user, so includes checks
#
# returns data.frame in requested parametrization
#
transform.df.parms = function(parms.df, out.format = c('b','beta','eta'), include.zero = TRUE)
{
  # start with many checks
  out.format = match.arg(out.format)
  colnames(parms.df) = tolower(colnames(parms.df))
  parms.df=ungroup(parms.df)
  
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
  
  parms.df = parms.df %>% 
    mutate(item_id = as.character(.data$item_id), item_score = as.integer(.data$item_score)) %>%
    arrange(.data$item_id, .data$item_score)
  
  
  mm = parms.df %>% 
    group_by(.data$item_id) %>% 
    summarise(min_score = min(.data$item_score), max_score = max(.data$item_score)) %>%
    ungroup()
  
  in.zero = any(mm$min_score == 0)
  
  if(in.zero && any(mm$min_score) != 0)
    stop("Either all items or none of the items should include a zero score parameter")
  
  if(any(mm$max_score == 0))
    stop('All items should contain at least one non-zero score parameter')
  
  if(any(mm$min_score<0))
    stop("Negative scores are not allowed")
  
  if(in.format == 'b' && any(parms.df$b < 0))
    stop("A 'b' parameter cannot be negative, perhaps you meant to include a 'beta' parameter?")
  
  fl = parms.df %>%
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>% 
    summarize(first = as.integer(min(.data$rn)), last = as.integer(max(.data$rn))) %>%
    ungroup() %>%
    arrange(.data$item_id)
  
  args = list(first = fl$first, last = fl$last, parms.df = parms.df, 
              out.zero = include.zero, in.zero = in.zero)
  do.call(get(paste0(in.format,'2',out.format)), args)[,c('item_id','item_score',out.format)] %>%
    arrange(.data$item_id)
}




#######################################################################
## Some functions to transform from one parameterization/Normalization to the other
#######################################################################

# If the zero category is present and has to be removed from
# parameters b, item_scores a, and index vectors first and last
remove_zero = function(a,b,first,last)
{
  if (is.matrix(b)) {
    b=b[,-first]
    if (is.null(dim(b))) b=as.matrix(t(b))
  }
  if (is.vector(b)) b=b[-first]
  a=a[-first]
  new_first=first
  new_last=last-1L
  for (i in 2:length(first))
  {
    ncat=last[i]-first[i]
    new_first[i]=new_last[i-1]+1L
    new_last[i]=new_first[i]+ncat-1L
  }
  return(list(a = a, b = b, first=new_first, last=new_last))
}

## add zero category 
add_zero = function(a, b, first,last)
{
  new_first=first[1]
  new_last=last[1]+1L
  if (length(first)>1)
  {
    for (i in 2:length(first))
    {
      nn=last[i]-first[i]
      new_first[i]=new_last[i-1]+1L
      new_last=c(new_last,new_first[i]+nn+1L)
    }
  }
  new_a = integer(length(a)+length(first))
  new_b = double(length(a)+length(first))
  new_a[-new_first]=a 
  new_b[new_first]=1
  new_b[-new_first]=b
  return(list(b=new_b, a=new_a, first=new_first, last=new_last))
}

## Makes the reparameterization matrix from log(b) to beta
makeD <- function(a,first,last)
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


################################################################
## Functions to go from one set of parameters to an other
# These are low-level functions with vectors or scalars as input
# and as output. Use apply when the input is a matrix.
#################################################################
beta2eta_ <-function(a, beta, first, last)
{
  eta = rep(0,length(beta))
  nI = length(first)
  for (i in 1:nI)
  {
    m = last[i]-first[i]
    eta[first[i]] = beta[first[i]]*a[first[i]]
    if (m>0)
    {
      for (j in (first[i]+1):last[i]) 
      {
        eta[j] = eta[j] + beta[first[i]]*a[first[i]]
        for (g in (first[i]+1):j)
        {
          eta[j] = eta[j] + beta[g]*(a[g]-a[g-1])
        }
      }
    }
  }
  return(eta)
}

eta2beta_ <-function(a, eta, first, last)
{
  beta = rep(0,length(eta))
  for (i in 1:length(first))
  {
    beta[first[i]] = eta[first[i]]/a[first[i]]
    for (j in (first[i]+1):last[i]) beta[j] = (eta[j]-eta[j-1])/(a[j]-a[j-1])
  }
  return(beta)
}

eta2b_ <- function(eta){exp(-eta)}

beta2b_ <-function(a,beta,first,last)
{
  eta = beta2eta_(a,beta,first,last)
  eta2b_(eta)
}

b2beta_ <-function(a,b,first,last)
{
  DD = makeD(a,first,last)
  beta = DD%*%log(b)
  return(beta)
}

b2eta_ <-function(a,b,first,last)
{
  DD = makeD(a,first,last)
  beta = DD%*%log(b)
  eta = beta2eta_(a, beta, first, last)
  return(eta)
}
####################################################################

### Change parameterization and normalization to produce OPCML output
## assumes that the first parameter is the reference unless there are fixed parameters
# TO DO: TB. Aanpassen zodat ook lambda wordt ge-renormaliseerd.
toOPLM = function(a, b, first, last, H=NULL, fixed_b=NULL, lambda=NULL)
{
  b_rn = b
  a_org = a
  if (!is.null(H)) H=H[-first,-first]
  if (!is.null(fixed_b)) fixed_b=fixed_b[-first]
  tmp = remove_zero(a,b,first,last)
  b = tmp$b; a = tmp$a
  first = tmp$first; last = tmp$last
  
  logb=log(b)
  cov.beta=NULL
  ########################
  
  ### Bayesian: b is a matrix
  if (is.matrix(b))
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
        c = mean(beta[r,])
        b_rn[r,] = b_rn[r,]*exp(c*a_org)
        beta[r,] = beta[r,] - c ## mean center
      }
    }
    if (nrow(beta)>2) cov.beta=cov(beta)
  }else
  {                                       ### CML; b is a vector
    DD = makeD(a,first,last)
    if (is.null(fixed_b))
    {
      beta = DD%*%logb
      b_rn = b_rn*exp(mean(beta)*a_org) # re-normalize b such that it corresponds to beta
      k  = length(b)
      CC = matrix(-1/k,k,k); diag(CC)=(k-1)/k
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
        fixed_set = which(!is.na(fixed_b))
        cov.beta  = solve(H[-fixed_set,-fixed_set])
        cov.beta  = DD[,-fixed_set]%*%cov.beta%*%t(DD[,-fixed_set])
      }
    }
  }  
  return(list(beta=beta, cov.beta=cov.beta, a=a, b_renorm = b_rn, first=first, last=last))
}

## Thus function expects category thresholds beta, a vector of item_category scores a,
#  and first and last. All without the zero category.
#  It returns dexter parameters b, as well as new a, first and last with the zero category.
toDexter <- function(beta, a, first, last, re_normalize=TRUE)
{
  if (re_normalize) beta = beta - beta[1]  # normalize different
  b = beta2b_(a,beta,first,last)
  names(b)=names(beta)
  
  ## add zero category
  tmp = add_zero(a,b,first,last)
  
  ## put everything in a (minimal) parms object
  est=list(b=tmp$b, a=tmp$a)
  inputs=list(ssIS=list(item_score=tmp$a),ssI=list(first=tmp$first,last=tmp$last))
  parms = list(est=est, inputs=inputs)
  return(parms)
}


#####################
# Some functions to transform user-provided (i.e., fixed) parameter values from one 
# parameterization to the other.
#####################
beta2eta <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$eta = beta2eta_(df.new$item_score, df.new$beta, first, last)
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(parms.df$item_score, parms.df$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          beta = rep(0,length(tmp$a)), eta = rep(0,length(tmp$a)),
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$beta[(tmp$first[i]+1):tmp$last[i]] = parms.df$beta[first[i]:last[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

beta2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$b = eta2b_(beta2eta_(df.new$item_score, df.new$beta, first, last))
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(df.new$item_score, df.new$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          beta = rep(0,length(tmp$a)), b = tmp$b,
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$beta[(tmp$first[i]+1):tmp$last[i]] = parms.df$beta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

eta2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$b = eta2b_(df.new$eta)
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(df.new$item_score, df.new$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          eta = rep(0,length(tmp$a)), b = tmp$b,
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

b2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=TRUE)
{
  df.new=parms.df
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(parms.df$item_score, parms.df$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          b = tmp$b, stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

