# common utility functions and datasets


# kind of fit a monotone B-spline
# a is an object containing Timo's EAPs, output of ability_tables
# monosm = function(a){
#   nss = nrow(a)
#   if (nrow(a)<4) return(a$theta) # no smoothing if too short
#   x = a$sumScore
#   y = a$theta
#   typical_no_knots = 15
#   typical_order = 5
#   nik = min(nss, typical_no_knots)
#   order = typical_order
#   innerknots = seq(x[1], x[nss], length=nik)
#   multiplicities = rep(1, nik-2)
#   lowend = x[1]
#   highend = x[nss]
#   innerknots = innerknots[2:(nik-1)]
#   knots = extendPartition (innerknots, multiplicities, order,
#                            lowend, highend)$knots
#   h = bsplineBasis (x, knots, order)
#   g = rowSums(h) - t(apply (h, 1, cumsum))
#   g = cbind (1, g)
#   u = pnnls (g, x, 1)$x
#   v = g%*%u
#   return(v)
# }

# differs in result from ntile in taking weights and in never putting equal values in different bins
# take care that if the nbr of distinct values is close to n, this will lead to very unequal sized bins
weighted_ntile = function(x, weights, n)
{
  
  rn = tibble(x=x,w=weights,ord=1L:length(x)) %>%
    arrange(.data$x) %>%
    mutate(rn=cumsum(.data$w)-.data$w) %>%
    arrange(.data$ord) %>%
    pull(.data$rn)
  
  len = sum(weights)
  as.integer(floor(n * rn/len + 1))
}

# other option, less memory efficient and does split equal values
# weighted_ntile2 = function(x, weights, n)
# {
#   x = rep(x,weights)
#   ntile(x,n)[cumsum(weights)-weights+1]
# }


# does basic argument type and attribute checks with error messages
# to do:
# one of multiple possible types

check_arg = function(x, type, name = deparse(substitute(x)), nullable = FALSE, .length = NA )
{
  if(is.null(x))
  {
    if(!nullable)
      stop(paste0("Argument'",name, "' may not be NULL"))
    
    return(NULL)
  }
  
  if(type == 'dataSrc')
  {
    if(!(inherits(x, 'dx_resp_data') || inherits(x, 'data.frame') || inherits(x, 'DBIConnection')))
    {
      stop(paste0("Argument'",name, "' must be of type 'DBIConnection' or 'data.frame'"))
    }
  } else if(type == 'integer')
  {
    if(!is.numeric(x) || x%%1 != 0)
      stop(paste0("Argument'",name, "' must be an integer value"))
    
  } else if(type %in% c('numeric','double'))
  {
    if(!is.numeric(x))
      stop(paste0("Argument'",name, "' must be numeric"))
    
  } else if(!inherits(x, type))
  {
      stop(paste0("Argument'",name, "' must be of type '", type, "'"))
  }
  
  if(!is.na(.length) && length(x) != .length)
    stop(paste0("Argument'",name, "' must have length ", .length))
}


dropNulls = function (x) 
{
  x[!vapply(x, is.null, FUN.VALUE = logical(1))]
}

# each time counter is read using $get, it increases by one
counter = setRefClass('counter',
  fields=list(x='numeric'), 
  methods=list(
    initialize = function(...)
    {
      callSuper(..., x=0)
    },
    get=function()
    { 
      x <<- x+1;
      return(x)
    }))

# use for forwarding arguments to e.g. plot function
merge_arglists = function(args, default = NULL, override = NULL)
{
  if(is.null(default))
    default = list()
  
  if(is.null(override))
    override = list()
  
  res = modifyList(default, args)
  res = modifyList(res, override)
  res
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

possible_scores <- function(a,first,last)
{
  y=a[first[1]:last[1]]
  for (i in 2:length(first))
  {
    y=sort(unique(as.vector((outer(y,a[first[i]:last[i]],'+')))))
  }
  return(y)
}


first_last2indx = function(first,last) unlist(apply(data.frame(first,last),1,function(x) x[1]:x[2]))

# internal utility function
# @parameter ssI: ssI as found in parms object
# @parameter design: data.frame with a column item_id
# we assume design is valid(i.e. doesn't contain items not existing in ssI)
# @return ssI like in parms but for a single booklet
subset_ssI = function(ssI, design)
{	
  # I see no need to recompute nCat, since it's already there in ssI
  ssI %>%
    semi_join(design, by='item_id') %>%
    arrange(.data$item_id) %>%
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1,last = cumsum(.data$nCat))	 
}


# highest posterior density interval
# not safe for bimodal distributions
hpd=function(x, conf=0.95, print=TRUE)
{
  conf <- min(conf, 1-conf)
  n <- length(x)
  nn <- round( n*conf )
  x <- sort(x)
  xx <- x[ (n-nn+1):n ] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx==m)[1]
  if (print)
  {
    return(data.frame(l=x[ nnn ],r=x[ n-nn+nnn ]))
  }else
  {
    return(rbind(x[ nnn ],x[ n-nn+nnn ]))
  }
}

# For a discrete vector x, this function gives 
# Frequencies P(x operator i) for i in min_max[1]:min_max[2]
# where operator can be ==, <, <=, etc.
# is min_max=NULL the range is min(x):max(x)
my_freq = function(x, min_max=NULL, operator=c("==", "<", "<=", ">", ">="))
{
  operator = match.arg(operator)
  if (!operator%in%c("==", "<", "<=", ">", ">=")) stop("wrong operator in my_freq")
  
  if (is.null(min_max)){
    rng = min(x):max(x)
  }else
  {
    rng = min_max[1]:min_max[2]
  }
  
  out = NULL
  if (operator=="==")
  {
    for (i in rng)
    {
      out = c(out,sum(x==i))
    }
  }
  
  if (operator=="<")
  {
    for (i in rng)
    {
      out = c(out,sum(x<i))
    }
  }
  
  if (operator=="<=")
  {
    for (i in rng)
    {
      out = c(out,sum(x<=i))
    }
  }
  
  if (operator==">")
  {
    for (i in rng)
    {
      out = c(out,sum(x>i))
    }
  }
  
  if (operator==">=")
  {
    for (i in rng)
    {
      out = c(out,sum(x>=i))
    }
  }
  
  out = out/length(x)
  names(out) = as.character(rng)
  return(out)
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
  out=gg$vectors[,which.max(gg$values)]
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
  kfun <- function(x, k1, k2) sort(x, partial = sort(c(k1, 
                                                       k2)))[c(k1, k2)]
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



#########################
#' Verbal aggression data
#' 
#' A data set of self-reported verbal behaviour in different frustrating
#' situations (Vansteelandt, 2000)
#' 
#' 
#' @name verbAggrData
#' @docType data
#' @format A data set with 316 rows and 26 columns.
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
#' @format A data set with 72 rows and 3 columns (item_id, response, item_score).
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





