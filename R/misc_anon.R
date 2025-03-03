


# @param setA, setB: two mutually exclusive subsets of items as indexes in first/last
# @return         a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b
# if cIM is not null, the interaction model will be used
# cIM should be per score
SSTable = function(b, a, first, last,setA, setB, cIM_score=NULL)
{
  firstA = first[setA] -1L
  firstB = first[setB] -1L
  lastA = last[setA] -1L
  lastB = last[setB] -1L
  
  if(!is.null(cIM_score))
  {
    sstable_imC(a, b, cIM_score, firstA, lastA, firstB, lastB)
  } else
  {
    sstable_nrmC(a, b, firstA, lastA, firstB, lastB)
  }
}




### Greatest Common Divisor via Euclid's algorithm
GCD2_ =function (n, m) 
{
  if (n == 0 && m == 0) 
    return(0)
  n = abs(n)
  m = abs(m)
  if (m > n) {
    t = n
    n = m
    m = t
  }
  while (m > 0) {
    t = n
    n = m
    m = t%%m
  }
  return(n)
}


GCD_ =function (x) 
{
  stopifnot(is.numeric(x))
  if (floor(x) != ceiling(x) || length(x) < 2) 
    stop("Argument 'x' must be an integer vector of length >= 2.")
  x = x[x != 0]
  n = length(x)
  if (n == 0) {
    g = 0
  }
  else if (n == 1) {
    g = x
  }
  else if (n == 2) {
    g = GCD2_(x[1], x[2])
  }
  else {
    g = GCD2_(x[1], x[2])
    for (i in 3:n) {
      g = GCD2_(g, x[i])
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
  conf = min(conf, 1-conf)
  n = length(x)
  nn = round( n*conf )
  x = sort(x)
  xx = x[ (n-nn+1):n ] - x[1:nn]
  m = min(xx)
  nnn = which(xx==m)[1]
  return(c(l=x[ nnn ],r=x[ n-nn+nnn ]))
}



# This function calculates overall and pointwise confidence envelopes 
# for a curve based on replicates of the curve evaluated at a number of fixed points.
# Based on theory by Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. 
# Cambridge University Press. Insprired by code from package boot.

# mat is a matrix with nrow = nr of replications, ncol = nr of points
# Example: test information for each of ncol ability values is calculated for nrow samples of 
# item parameters from posterior.

#TO~DO: protect against NA's
conf_env = function(mat, level = 0.95) 
{
  overall_found = TRUE
  emperr = function(rmat, p = 0.05, k = NULL) {
    R = nrow(rmat)
    if (is.null(k)) 
      k = p * (R + 1)/2
    else p = 2 * k/(R + 1)
    kf = function(x, k, R) 1 * ((min(x) <= k) | (max(x) >= 
                                                    R + 1L - k))
    c(k, p, sum(apply(rmat, 1L, kf, k, R))/(R + 1))
  }
  kfun = function(x, k1, k2) sort(x, partial = sort(c(k1, k2)))[c(k1, k2)]
  index = 1L:ncol(mat)
  if (length(index) < 2L) 
    stop("This function for curves")
  rmat = apply(mat, 2L, rank)
  R = nrow(mat)
  if (length(level) == 1L) 
    level = rep(level, 2L)
  k.pt = floor((R + 1) * (1 - level[1L])/2 + 1e-10)
  k.pt = c(k.pt, R + 1 - k.pt)
  err.pt = emperr(rmat, k = k.pt[1L])
  ov = emperr(rmat, k = 1)
  ee = err.pt
  al = 1 - level[2L]
  if (ov[3L] > al) 
    overall_found = FALSE
  else {
    continue = !(ee[3L] < al)
    while (continue) {
      kk = ov[1L] + round((ee[1L] - ov[1L]) * (al - ov[3L])/(ee[3L] - 
                                                                ov[3L]))
      if (kk == ov[1L]) 
        kk = kk + 1
      else if (kk == ee[1L]) 
        kk = kk - 1
      temp = emperr(rmat, k = kk)
      if (temp[3L] > al) 
        ee = temp
      else ov = temp
      continue = !(ee[1L] == ov[1L] + 1)
    }
  }
  k.ov = c(ov[1L], R + 1 - ov[1L])
  err.ov = ov[-1L]
  out = apply(mat, 2L, kfun, k.pt, k.ov)
  if (overall_found){
    out = out[4:3, ]
  }else
  {
    out = out[1:2, ]
  }
  return(out)
}


# equivalent to log(sum(exp(x)))
# where sum(exp(x)) is potentially infinite in floating point
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
c2weights = function(cIM)
{
  hh=kronecker(t(cIM),cIM,'+')
  gg=eigen(hh)
  out=abs(gg$vectors[,which.max(gg$values)])
  av_indx=which.min(abs(out-mean(out)))
  return(out/out[av_indx])
}


all_trivial_scores = function(scores)
{
  s = lapply(split(scores$item_score, scores$item_id, drop=TRUE),function(x) c(0L,x))
  length(Reduce(function(a,b){out = as.vector(outer(a,b,'+')); if(!anyDuplicated(out)) out else NULL},s)) > 0
}


# differs in result from ntile in taking weights and in never putting equal values in different bins
# integer weights only
weighted_ntile = function(x, weights, nbins)
{
  nbins = as.integer(nbins)
  weights = as.integer(weights)
  
  if(length(x) <= nbins)
    return(as.integer(rank(x)))
  
  if(is.unsorted(x))
  {
    dat = tibble(x=x, w=weights, ord=1:length(x)) |>
      arrange(.data$x)
    
    dat$bin = weighted_binning(dat$w,nbins)
    dat |> arrange(ord) |> pull(bin) |> drop()
  } else
  {
    drop(weighted_binning(weights,nbins))
  }
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





