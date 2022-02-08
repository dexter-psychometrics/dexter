# common utility functions


.onLoad = function(libname, pkgname) 
{
  dexter_options = list(dexter.use_tibble=FALSE, dexter.progress=TRUE)
  
  to_set = setdiff(names(dexter_options), names(options()))
  if(length(to_set)>0)
    options(dexter_options[to_set])

  invisible()
}


explicit_NA = function(x, replace_NA_with = c('<NA>','.NA.','<<NA>>','__NA__'))
{
  if(!is.character(x) || !anyNA(x))
    return(x)
  
  for(v in replace_NA_with)
  {
    if(!v %in% x)
    {
      x[is.na(x)] = v
      return(x)
    }
  }
  stop('could not resolve NA values')
  
}


## Burnin and thinning for different purposes:
# cal: Bayesian calibration
# pv: plausible values
# ps: plausible scores
Gibbs.settings = list(from.cal = 20L, step.cal = 2L, start_b='cml',
                      from.pv = 20L, step.pv = 5L,
                      from.ps = 1L, step.ps = 1L)



which_gibbs = function(n_samples, n_gibbs)
{
  step = n_gibbs %/% n_samples
  from = (n_gibbs %% n_samples) + 1L
  
  which = 
    if(step>0)
      seq(from, n_gibbs, step)
    else
      1:n_samples
  
  list(step=step, from=from, which=which, enough = step>0)
}


df_format = function(df)
{
  if(getOption('dexter.use_tibble', FALSE))
  {
    as_tibble(ungroup(df))
  } else
  {
    as.data.frame(df, stringsAsFactors = FALSE)
  }
}



is.date = function(x) inherits(x, "Date")
is.time = function(x) inherits(x,'POSIXt')

# add named column(s) to the end of a data.frame
add_column = function(df, ...)
{
  dots = list(...)
  if(is.null(names(dots)) || any(names(dots) == ''))
    stop('arguments must be named')
  
  for(nm in names(dots))
  {
    vec = dots[[nm]]
    if(NROW(df)==0)
      vec = vec[0]
    
    stopifnot(length(vec)==1 || NROW(df) == length(vec))
    df[[nm]] = vec
  }
  
  df
}

bind_vertical = function(lst)
{
  mtx = any(sapply(lst, is.matrix))
  if(mtx)
    bind_rows(lst)
  else
    matrix(unlist(lst), ncol=1)
}


# scores is data.frame item_id, item_score, assumed that 0 is present
all_trivial_scores = function(scores)
{
  length(Reduce(function(a,b){out = as.vector(outer(a,b,'+')); if(!anyDuplicated(out)) out else NULL}, 
                split(scores$item_score, scores$item_id, drop=TRUE))
  ) > 0
}




# format string with named arguments

# txt: format string with named arguments prefixed by a dollar sign, formatting can be done with postfixing with : 
# arglist: list with named arguments to interpolate in the format string. Use only alphanumerical characters in the names

# return: 
# formatted string

# examples

# fstr('$bla, $b',list(bla='some string'))
# fstr('$bla:.1f, $b',list(bla=3.2423))
#
fstr = function(txt, arglist)
{
  if(length(txt) != 1 || length(arglist)==0 || is.null(txt) || !grepl('$',txt,fixed=TRUE)) 
    return(txt)

  txt = gsub('%','%%',txt,fixed=TRUE)
  
  spr_args = list()
  txt_copy = txt
  matches = gregexpr('\\$[a-z]\\w*(\\:[\\.\\d]*[ifscega])?', txt, ignore.case = TRUE, perl=TRUE)[[1]]
  matches_end = attr(matches, 'match.length') + matches - 1L
  for(i in seq_along(matches))
  {
    matched_str = substring(txt, matches[i], matches_end[i])
    matched_tpl = unlist(strsplit(matched_str,':',fixed=TRUE))
    matched_name = substring(matched_tpl[1],2)
    # text can come after a match if match not ends in :
    if(length(matched_tpl) == 1 && !(matched_name %in% names(arglist)))
    {
      l = -1
      m = NULL
      for(nm in names(arglist))
      {
        if(startsWith(matched_name, nm) && nchar(nm) > l)
        {
          l = nchar(nm)
          m = nm
        }
      }
      if(!is.null(m))
      {
        matched_name = m
        matched_str = paste0('$',m)
      }
    }
    
    if(matched_name %in% names(arglist))
    {
      arg = arglist[[matched_name]]
      if(length(matched_tpl) == 2)
      {
        txt_copy = sub(matched_str, paste0('%',matched_tpl[2]) ,txt_copy, fixed=TRUE)
      } else
      {
        arg = as.character(arg)
        txt_copy = sub(matched_str, '%s' ,txt_copy, fixed=TRUE)
      }
      spr_args = append(spr_args, arg) 
    }
  }
  
  do.call(sprintf, append(txt_copy, spr_args))
}


# differs in result from ntile in taking weights and in never putting equal values in different bins
# take care that if the nbr of distinct values is close to n, this will lead to very unequal sized bins
weighted_ntile = function(x, weights, n)
{
  
  dat = tibble(x=x, w=weights, ord=1:length(x)) %>%
    arrange(.data$x) %>%
    mutate(rn=cumsum(.data$w)-.data$w) %>%
    arrange(.data$ord) 
  
  as.integer(floor(n * dat$rn/sum(dat$w) + 1))
}

weighted_cor = function(x,y,n)
{
  x_u = x - weighted.mean(x,n)
  y_u = y - weighted.mean(y,n)

  weighted.mean(x_u * y_u, n)/
    (sqrt(weighted.mean(x_u ^ 2, n)) * sqrt(weighted.mean(y_u ^ 2, n)))
}

# not able to deal with na,nan,etc
weighted_quantile = function(x,w,probs)
{
  if(is.unsorted(x))
  {
    ord = order(x)
    x = x[ord]
    w = w[ord]
  }  
  csw = cumsum(w)
  n = sum(w)
  #np = length(probs)
  index = 1 + (n - 1) * probs
  lo = floor(index)
  hi = ceiling(index)
  
  qs = x[sapply(lo,function(i) min(which(csw>=i)))]
  
  i = which(index > lo)
  h = (index - lo)[i]
  qs[i] = (1 - h) * qs[i] + h * x[sapply(hi,function(i) min(which(csw>=i)))]
  
  qs
  
}

# variance of total sample by combining group variances
combined_var = function(means,vars,n)
{
  if(length(vars)<=1L)
    return(vars)
  q = (n-1)*vars + n*means^2
  (sum(q) - sum(n)* weighted.mean(means,n)^2)/(sum(n)-1)
}



# non vectorized version of ifelse
if.else = function(test, yes, no)
{
  if(isTRUE(test)) return(yes)
  no
}


#  basic argument type and attribute checks with error messages
stop_ = function(...) stop(..., call. = FALSE)

check_dataSrc = function(x)
{
  force(x)
  if(inherits(x, 'dx_resp_data'))
     return(NULL)
  
  if(is.matrix(x))
  {
    if(!mode(x) %in% c('numeric','integer'))
      stop_('dataSrc must be a matrix of positive numbers, a data.frame or a database connection')
    return(NULL)
  } 
  
  if(inherits(x, 'data.frame'))
  {
    if(length(setdiff(c('person_id','item_id','item_score'), colnames(x)))>0)
      stop_("dataSrc must contain the columns: 'person_id','item_id','item_score'")
    return(NULL)
  }
  if(is_db(x))
  {
    if(dbIsValid(x)) return(NULL)
    stop_('your database connection is no longer valid, you need to reconnect. see: ?open_project for details')
  }
    
     
  if(length(x)== 1 && is.character(x) && file.exists(x))
      stop_('dataSrc is a string but must be a database connection, data.frame or matrix. ',
           'Did you forget to do: `db = open_project("',x,'")`?')
 
  stop_("dataSrc must be of type 'DBIconnection', 'data.frame' or 'matrix'")
}

check_db = function(x)
{
  if(length(x)== 1 && is.character(x) && file.exists(x))
  {
      stop_('db is a string but must be a database connection. ',
            'Did you forget to do: `db = open_project("',x,'")`?')
  } else if(!is_db(x))
  {
    stop_("Argument 'db' is not a database connection.")
  } else if(!dbIsValid(x))
  {
    stop_('your database connection is no longer valid, you need to reconnect. see: ?open_project for details')
  }
  
  NULL
}


check_vector = function(x, type=c('character','numeric','integer'), name = deparse(substitute(x)), 
                        nullable = FALSE, .length = NULL, .min=NULL, .max=NULL )
{
  if(nullable && is.null(x))
    return(NULL)
  
  type = match.arg(type)
  
  if(!is.vector(x) || is.list(x))
    stop("Argument '",name, "' must be a vector of type ", type, call.=FALSE)
  
  if(type != 'character')
  {
    if(!is.numeric(x))
      stop_("Argument '",name, "' must be ", type)
  
    if(type=='integer' && !(is.integer(x) || all(x %% 1 == 0)))
       stop_("Argument '",name, "' must be an integer")
    
    if(!is.null(.min) && any(x<.min))
      stop_("Argument '",name, "' must be >= ", .min)
    
    if(!is.null(.max) && any(x>.max))
      stop_("Argument '",name, "' must be <= ", .max)
  }
  

  if(!is.null(.length) && length(x) != .length )  
    stop_("Argument '",name, "' must have length ", .length)
 
}
  
  
check_num = function(x, type=c('numeric','integer'), name = deparse(substitute(x)), 
                     nullable = FALSE, .length = NULL, .min=NULL, .max=NULL )
{
  type=match.arg(type)
  name=force(name)
  check_vector(x,type=type,name=name,nullable=nullable,.length=.length,.min=.min,.max=.max)
}


check_character = function(x, name = deparse(substitute(x)), nullable = FALSE, .length = NULL)
{
  name=force(name)
  check_vector(x,type='character',name=name,nullable=nullable,.length=.length)
}

check_string = function(x, name = deparse(substitute(x)), nullable = FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  
  if(!is.character(x) && length(x)!=1)
    stop_("Argument '",name, "' must be a string of length 1")
}


check_df = function(x, columns=NULL, n_rows=NULL, name = deparse(substitute(x)), nullable=FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  
  if(!inherits(x, 'data.frame'))
    stop("Argument'",name, "' must be a data.frame")
  
  missing_col = setdiff(columns, colnames(x))
  
  if(!is.null(columns) && length(missing_col>0))
    stop('column(s): ', paste0('`', missing_col, '`',collapse=', '),' must be present in ', name)
  
  if(!is.null(n_rows) && NROW(x)!=n_rows)
    stop_('argument`', name, '` must have ', n_rows, ' rows')
  
}

check_parms = function(x, name = deparse(substitute(x)), nullable=FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  if(!inherits(x,'prms'))
    stop_(name,' must be an object of type `prms`')
}

check_list = function(x, name = deparse(substitute(x)), nullable=FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  if(!inherits(x,'list'))
    stop_(name,' must be a list')
}


# start, stop
# rng_fl(1:6) => c(1,6)
# rng_fl(c(1,6)) => c(1,6)
# rng_fl(c(6,1)) => error, etc.
rng_fl = function(x, name = deparse(substitute(x)))
{
  if(length(x)==2)
  {
    if(x[1]>x[2])
      stop_('first element of ',name,' must be smaller or equal than second element')
    x
  } else if(length(x)<2)
  {
    stop_(name, ' must be a vector of length 2')
  } else
  {
    test = x[1]:x[length(x)] 
    if(length(x) != length(test) || !all(test == x) || x[1]>x[length(x)])
      stop_(name, ' is not a valid range')
    c(x[1],x[length(x)])
  }
}




dropNulls = function (x) 
{
  x[!vapply(x, is.null, FUN.VALUE = logical(1))]
}

# use for forwarding arguments to e.g. plot function
merge_arglists = function(args, default = NULL, override = NULL)
{
  if(!is.null(default))
    args = modifyList(default, args)
  
  if(!is.null(override))
    args = modifyList(args, override)

  args
}


df_identical = function(a, b)
{
  # check all values in dataframe equal, disregard column order
  
  if(!all(dim(a)==dim(b))) return(FALSE)
  if(!length(intersect(colnames(a),colnames(b))) == ncol(a)) return(FALSE)
  
  a = a %>% mutate_if(is.factor, as.character) 
  b = b %>% mutate_if(is.factor, as.character)
  
  return(all(a == b[,colnames(a)]))
}

is_connected = function(design)
{
  design = droplevels(design)
  # usually best via items but with extreme predicates and large data this might lead to 
  # adj larger than working memory
  if(nlevels(design$item_id) >= nlevels(design$booklet_id))
  {
    items = as.matrix(table(design$item_id, design$booklet_id))
    adj = crossprod(items, items)
  } else
  {
    booklets = as.matrix(table(design$booklet_id, design$item_id))
    adj = crossprod(booklets, booklets)
  }
    
  mode(adj) = 'integer'
  
  all(ds_connected_groups(adj)==1)
}

### Greatest Common Divisor via Euclid's algorithm
GCD2_ <-function (n, m) 
{
  if (n == 0 && m == 0) 
    return(0)
  n <- abs(n)
  m <- abs(m)
  if (m > n) {
    t <- n
    n <- m
    m <- t
  }
  while (m > 0) {
    t <- n
    n <- m
    m <- t%%m
  }
  return(n)
}


GCD_ <-function (x) 
{
  stopifnot(is.numeric(x))
  if (floor(x) != ceiling(x) || length(x) < 2) 
    stop("Argument 'x' must be an integer vector of length >= 2.")
  x <- x[x != 0]
  n <- length(x)
  if (n == 0) {
    g <- 0
  }
  else if (n == 1) {
    g <- x
  }
  else if (n == 2) {
    g <- GCD2_(x[1], x[2])
  }
  else {
    g <- GCD2_(x[1], x[2])
    for (i in 3:n) {
      g <- GCD2_(g, x[i])
      if (g == 1) 
        break
    }
  }
  return(g)
}

# highest posterior density interval
# not safe for bimodal distributions
hpdens = function(x, conf=0.95)
{
  conf <- min(conf, 1-conf)
  n <- length(x)
  nn <- round( n*conf )
  x <- sort(x)
  xx <- x[ (n-nn+1):n ] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx==m)[1]
  return(c(l=x[ nnn ],r=x[ n-nn+nnn ]))
}


# For a discrete vector x, this function gives 
# Frequencies P(x operator i) for i in min_max[1]:min_max[2]
# where operator can be ==, <, <=, etc.
# is min_max=NULL the range is min(x):max(x)
my_freq = function(x, min_max = min(x):max(x), operator = c("==", "<", "<=", ">", ">=","!="))
{
  if(!is.function(operator))
    operator = get(match.arg(operator))
  
  out = sapply(min_max, function(i) {sum(operator(x, i))} )/length(x)
  names(out) = min_max
  out
}

geo_mean = function(x)
{
   return(exp(mean(log(x))))
}
### round to geometric mean ###
# Round positive numbers x to integers
# such that the geometric mean equals approximately J
r2gm = function(x, J)
{
  G=geo_mean(x)
  out=NULL
  for (i in 1:length(x)) out = c(out,max(1,floor(0.5+(J*x[i]/G))))
  return(out)
}

### Weights based on a rank-one approximation to the
# Interaction model
c2weights<-function(cIM)
{
  hh=kronecker(t(cIM),cIM,'+')
  gg=eigen(hh)
  out=abs(gg$vectors[,which.max(gg$values)])
  av_indx=which.min(abs(out-mean(out)))
  return(out/out[av_indx])
}

#########
# This function calculates overall and pointwise confidence envelopes 
# for a curve based on replicates of the curve evaluated at a number of fixed points.
# Based on theory by Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. 
# Cambridge University Press. Insprired by code from package boot.
##
# mat is a matrix with nrow = nr of replications, ncol = nr of points
# Example: test information for each of ncol ability values is calculated for nrow samples of 
# item parameters from posterior.
########
#TO~DO: protect against NA's
conf_env = function(mat, level = 0.95) 
{
  overall_found = TRUE
  emperr <- function(rmat, p = 0.05, k = NULL) {
    R <- nrow(rmat)
    if (is.null(k)) 
      k <- p * (R + 1)/2
    else p <- 2 * k/(R + 1)
    kf <- function(x, k, R) 1 * ((min(x) <= k) | (max(x) >= 
                                                    R + 1L - k))
    c(k, p, sum(apply(rmat, 1L, kf, k, R))/(R + 1))
  }
  kfun <- function(x, k1, k2) sort(x, partial = sort(c(k1, k2)))[c(k1, k2)]
  index = 1L:ncol(mat)
  if (length(index) < 2L) 
    stop("This function for curves")
  rmat <- apply(mat, 2L, rank)
  R <- nrow(mat)
  if (length(level) == 1L) 
    level <- rep(level, 2L)
  k.pt <- floor((R + 1) * (1 - level[1L])/2 + 1e-10)
  k.pt <- c(k.pt, R + 1 - k.pt)
  err.pt <- emperr(rmat, k = k.pt[1L])
  ov <- emperr(rmat, k = 1)
  ee <- err.pt
  al <- 1 - level[2L]
  if (ov[3L] > al) 
    overall_found = FALSE
  else {
    continue <- !(ee[3L] < al)
    while (continue) {
      kk <- ov[1L] + round((ee[1L] - ov[1L]) * (al - ov[3L])/(ee[3L] - 
                                                                ov[3L]))
      if (kk == ov[1L]) 
        kk <- kk + 1
      else if (kk == ee[1L]) 
        kk <- kk - 1
      temp <- emperr(rmat, k = kk)
      if (temp[3L] > al) 
        ee <- temp
      else ov <- temp
      continue <- !(ee[1L] == ov[1L] + 1)
    }
  }
  k.ov <- c(ov[1L], R + 1 - ov[1L])
  err.ov <- ov[-1L]
  out <- apply(mat, 2L, kfun, k.pt, k.ov)
  if (overall_found){
    out = out[4:3, ]
  }else
  {
    out = out[1:2, ]
  }
  return(out)
}

## Regular Block-Diagnal Matrix from List of matrices
# Courtesy of C.Ladroue
blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

#
# mean and variance of Y|x if x,y is multivariate normal
#
# with mean mu and variance-covariance matrix sigma
# @param m vector of means
# @param sigma covariance matrix
# @param y.ind indices dependent variable(s)
# @param x.ind indices conditioning variables. If null its just all others
# @param x.value value of conditioning variables
# 
condMoments = function(mu, sigma, y.ind, x.ind=NULL, x.value )
{
  if (is.null(x.ind)) x.ind = setdiff(1:length(mu), y.ind)
  B = sigma[y.ind, y.ind]
  C = sigma[y.ind, x.ind, drop = FALSE]
  D = sigma[x.ind, x.ind]
  CDinv = C %*% solve(D)
  if (is.vector(x.value))
  {
    cMu = c(mu[y.ind] + CDinv %*% (x.value - mu[x.ind]))
  }else
  {
    nP = nrow(x.value)
    cMu = rep(0,nP)
    for (i in 1:nP) cMu[i] = mu[y.ind] + CDinv %*% (x.value[i,] - mu[x.ind])
  }
  cVar = B - CDinv %*% t(C)
  return(list(mu=cMu, sigma=cVar))
}


# log(sum(exp(x)))
# where exp(x) is potentially infinite in floating point
logsumexp = function(x)
{
  m = max(x)
  m + log(sum(exp(x-m)))
}



# GH points
# library(statmod)
# GH = gauss.quad.prob(160,'normal',mu=0,sigma=1)
# quadpoints = tibble(nodes=GH$nodes,weights=GH$weights) |>
#   filter(weights>1e-60) |>
#   arrange(weights) |>
#   as.list()
# usethis::use_data(quadpoints, internal = TRUE)

