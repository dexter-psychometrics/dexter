# common utility functions and datasets

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



