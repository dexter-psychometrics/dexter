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



#########################
#' Verbal aggression data
#' 
#' A data set of self-reported verbal behaviour in different frustrating
#' situations (Vansteelandt, 2000)
#' 
#' 
#' @name verbAggrData
#' @docType data
#' @format A data set with 316 rows and 25 columns.
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

###########################################
#' Item properties in the PISA 2012 example
#' 
#' A data set of item properties in the PISA 2012 example (see the
#' help screen for function add_booklet)
#' 
#' 
#' @name PISA_item_class
#' @docType data
#' @format A data set with 109 rows and 6 columns.
#' @keywords datasets
NULL



