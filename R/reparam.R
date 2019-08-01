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
        eta[j] = eta[j] + beta[j-1]*a[j-1]
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
toOPLM = function(a, b, first, last, H=NULL, fixed_b=NULL)
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

