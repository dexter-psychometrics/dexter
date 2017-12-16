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



# logger for debugging
# add message with debug.log$send()
# retreive and clear all messages with debug.log$retrieve(), returns list

logger = setRefClass('logger', 
  fields = list(stack='list', msg='logical'), 
  methods = list(
    send = function(obj, name = 'no_name')
      {
        if(msg) message(obj)
        stack[[length(stack)+1]] <<-  obj
        names(stack)[length(stack)] <<- name 
        invisible(obj)
      },
    retrieve = function()
      {
        res = stack
        stack <<- list()
        return(res)
      }
  )
)
# global object
debug.log = logger$new(msg=FALSE)
           
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
  if(!is.null(default))
    for(nm in names(default)) 
      if(! nm %in% names(args)) 
        args[[nm]] = default[[nm]]
      
      if(!is.null(override))
        for(nm in names(override)) 
          args[[nm]] = override[[nm]]
        
        return(args)
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
# not save for bimodal distributions
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



