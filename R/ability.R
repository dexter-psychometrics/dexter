


#' Estimate abilities
#'
#' Computes estimates of ability for persons or for booklet scores
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and, 
#' depending on parametrization, a column named either beta/delta, eta or b
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param design A data.frame with columns item_id and optionally booklet_id. If the column booklet_id is not included, the score 
#' transformation table will be based on all items found in the design. If design is NULL
#' and parms is an enorm fit object the score transformation table will be computed based on the test design 
#' that was used to fit the items. 
#' @param method   Maximum Likelihood (MLE), Expected A posteriori (EAP) or Weighted Likelihood (WLE)
#' @param prior    If an EAP estimate is produced one can choose a normal prior or
#'                 Jeffreys prior; i.e., a prior proportional to the square root of test information.
#' @param parms_draw When parms is Bayesian, parms_draw can be the index of the posterior sample of the item 
#' parameters that will be used for generating abilities. If parms_draw='average', the posterior mean is used. 
#' @param mu Mean of the normal prior
#' @param sigma Standard deviation of the normal prior
#' @param merge_within_persons for persons who were administered multiple booklets, 
#' whether to provide just one ability value (TRUE) or one per booklet(FALSE)
#' 
#' @return 
#' \describe{
#'   \item{ability}{a data.frame with columns: booklet_id, person_id, booklet_score, theta and optionally se (standard error) }
#'   \item{ability_tables}{a data.frame with columns: booklet_id, booklet_score, theta and optionally se (standard error)}
#' }
#' 
#' @details MLE estimates of ability will produce -Inf and Inf estimates for
#' the minimum (=0) and the maximum score on a booklet. If this is undesirable, 
#' we advise to use WLE. The WLE was proposed by Warm (1989) to reduce bias in the MLE and is also known
#' as the Warm estimator.
#'
#' @examples
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules, ":memory:")
#' add_booklet(db, verbAggrData, "agg")
#' 
#' f = fit_enorm(db)
#' 
#' mle = ability_tables(f, method="MLE")
#' eap = ability_tables(f, method="EAP", mu=0, sigma=1)
#' wle = ability_tables(f, method="WLE")
#' 
#' plot(wle$booklet_score, wle$theta, xlab="test-score", ylab="ability est.", pch=19)
#' points(mle$booklet_score, mle$theta, col="red", pch=19,)
#' points(eap$booklet_score, eap$theta, col="blue", pch=19)
#' legend("topleft", legend = c("WLE", "MLE", "EAP N(0,1)"), 
#'         col = c("black", "red", "blue"), bty = "n",pch = 19)
#' 
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory. 
#' Psychometrika, 54(3), 427-450. 
#' 
ability = function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP","WLE"), prior=c("normal", "Jeffreys"), 
                   parms_draw='average', mu=0, sigma=4, merge_within_persons=FALSE)
{
  check_dataSrc(dataSrc)

  method = match.arg(method)
  prior = match.arg(prior) 
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  if(inherits(parms,'prms'))
  {
    parms_check = parms$inputs$ssIS[,c('item_id','item_score')]
  } else if(inherits(parms,'data.frame'))
  {
    parms_check = distinct(ungroup(parms), .data$item_id,.data$item_score)
  }

  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=TRUE, env=env, 
                           parms_check=parms_check, merge_within_persons = merge_within_persons)
  

  abl = ability_tables(parms=parms, design = respData$design, method = method, prior=prior, parms_draw = parms_draw, 
                       mu=mu, sigma=sigma)
  
  abl$booklet_id = ffactor(abl$booklet_id, levels = levels(respData$design$booklet_id))
  
  respData$x |> 
    inner_join(abl, by = c("booklet_id", "booklet_score")) |> 
    select(any_of(c('booklet_id', 'person_id', 'booklet_score', 'theta', 'se'))) |>
    mutate_if(is.factor, as.character) |>
    df_format()
}



#' @rdname ability
ability_tables = function(parms, design = NULL, method = c("MLE","EAP","WLE"), prior=c("normal", "Jeffreys"), 
                          parms_draw = 'average', mu=0, sigma=4)
{
  method = match.arg(method)
  prior = match.arg(prior) 
  if(is.numeric(parms_draw)) check_num(parms_draw,.length=1)
  else parms_draw = match.arg(parms_draw)
  
  if(method=='EAP' && prior=="normal")
  {
    check_num(mu, .length=1)
    check_num(sigma, .length=1, .min=0)
  }
  
  if (method=="EAP" && prior=="Jeffreys") method="jEAP"
  
  simple_parms = simplify_parms(parms, design, parms_draw) 
  b = simple_parms$b
  a = simple_parms$a
  
  estimate = switch(method, 
                    'MLE'  = function(.){ theta_MLE(b, a, .$first, .$last, se=TRUE) }, 
                    'EAP'  = function(.){ theta_EAP_GH(b, a, .$first, .$last, mu=mu, sigma=sigma) },
                    'jEAP' = function(.){ theta_jEAP(b, a, .$first, .$last, se=TRUE) },
                    'WLE' = function(.){ theta_WLE(b, a, .$first, .$last, se=TRUE) })
  
  
  # under the assumption that we always get theta's for the vector 0:max_test_score 
  simple_parms$design |> 
    group_by(.data$booklet_id) |>
    do({
      est = estimate(.)
      tibble(booklet_score=0:(length(est$theta)-1), theta = est$theta,se=est$se)
    }) |>
    ungroup() |>
    mutate_if(is.factor, as.character) |>
    df_format()
}


# Computes likelihood and test information for internal use
#
# For a vector of thetas it returns:
# l = a matrix (nbr of response cats * length of theta) of the likelihood or log-likelihood if log=TRUE
# I = a vector of the information function computed at each theta = sum(P'^2/PQ)
# J = a vector of the derivative of the information function at each theta
#
# The vector theta can be a set of quadrature points or estimates to compute their SE
#
# Note: can not deal with Inf or NA values in theta
IJ_ = function(b, a, first, last, theta, log=FALSE)
{
  nI = length(first)
  nT = length(theta)
  I = matrix(0, nT, nI)
  J = matrix(0, nT, nI)
  logFi = double(nT)
  
  a = as.integer(a)
  
  IJ_c(theta, b, a, as.integer(first-1L), as.integer(last-1L), I, J,logFi)
  
  scores = 0:sum(a[last])
  
  l = sweep(outer(scores,theta), 2, logFi, '-')
  if (!log) l = exp(l)
  return(list(I=rowSums(I), J=rowSums(J), l=l))
}



theta_MLE = function(b,a,first,last, se=FALSE)
{
  a = as.integer(a)
  theta = theta_mle_sec(b, a, as.integer(first-1L), as.integer(last-1L))[,1,drop=TRUE]
  
  sem = NULL
  if (se)
  {
    # use r indexed first last for IJ
    f = IJ_(b,a,first,last, theta)
    sem = c(NA, 1/sqrt(f$I), NA)
  }
  
  return(list(theta = c(-Inf,theta,Inf), se=sem))
}

theta_WLE = function(b,a,first,last, se=FALSE)
{
  a = as.integer(a)
  theta = theta_wle_sec(b, a, as.integer(first-1L), as.integer(last-1L))[,1,drop=TRUE]
  
  sem = NULL
  if (se)
  {
    # use r indexed first last for IJ
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  
  return(list(theta = theta, se=sem))
}




## EAP using Jeffrey's prior: aka propto sqrt(information)
# Uses a weighted average to integrate over a grid defined by:
# grid_from, grid_to and grid_length.
theta_jEAP = function(b, a, first,last, se=FALSE, grid_from=-6, grid_to=6, grid_length=101) 
{
  theta_grid = seq(grid_from, grid_to, length=grid_length)
  f = IJ_(b,a,first,last,theta_grid)
  prior=sqrt(f$I)
  w = sweep(f$l, 2, prior, '*')
  theta = apply(w, 1, function(x) weighted.mean(theta_grid, w=x))
  sem=rep(NA,length(theta))
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  return(list(theta=theta,se=sem))
}

# Expected distribution given a vector theta
# return matrix, ncol=length(theta), nrow=nscores
pscore = function(theta, b, a, first, last)
{
  g = elsymC(b, a, first-1L, last-1L)
  score = 0:(length(g)-1)
  p = sapply(theta, function(tht) log(g) + score*tht)
  
  exp(sweep(p,2,apply(p,2,logsumexp),`-`))
}

# Expected scores given one or more ability values theta
E_score = function(theta,b,a,first,last)
{
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  a = as.integer(a)
  
  Escore_C(theta, b, a, first, last)[,1,drop=TRUE]
}

rscore_item = function(theta,b,a,first,last)
{
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  a = as.integer(a)
  sampleNRM_itemC(theta, b, a, first, last)
}


# se is always returned (arg se is ignored)
theta_EAP_GH = function(b, a, first,last, se=TRUE, mu=0, sigma=4)
{
  nodes = quadpoints$nodes * sigma + mu
  weights = quadpoints$weights
  ps = t(pscore(nodes,b,a,first,last))
  lapply(theta_EAP_GH_c(ps,nodes,weights), drop)
}



#' Functions of theta
#' 
#' returns information function, expected score function, score simulation function, or score distribution 
#' for a single item, an arbitrary group of items or all items
#' 
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and, 
#' depending on parametrization, a column named either beta/delta, eta or b
#' @param items vector of one or more item_id's. If NULL and booklet_id is also NULL, all items in parms are used
#' @param booklet_id id of a single booklet (e.g. the test information function), if items is not NULL this is ignored
#' @param which.draw the number of the random draw (only applicable if calibration method was Bayes). If NULL, the mean 
#' beta parameter will be used
#' 
#' @return Each function returns a new function which accepts a vector of theta's. These return the following values: 
#' \describe{
#' \item{information}{an equal length vector with the information estimate at each value of theta.}
#' \item{expected_score}{an equal length vector with the expected score at each value of theta}
#' \item{r_score}{a matrix with length(theta) rows and one column for each item containing simulated scores based on theta. 
#' To obtain test scores, use rowSums on this matrix}
#' \item{p_score}{a matrix with length(theta) rows and one column for each possible sumscore containing the probability of 
#' the score given theta}
#' }
#' 
#' @examples
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules,':memory:')
#' add_booklet(db,verbAggrData, "agg")
#' p = fit_enorm(db)
#' 
#' # plot information function for single item
#' 
#' ifun = information(p, "S1DoScold")
#' 
#' plot(ifun,from=-4,to=4)
#' 
#' # compare test information function to the population ability distribution
#' 
#' ifun = information(p, booklet="agg")
#' 
#' pv = plausible_values(db,p)
#' 
#' op = par(no.readonly=TRUE)
#' par(mar = c(5,4,2,4))
#' 
#' plot(ifun,from=-4,to=4, xlab='theta', ylab='test information')
#' 
#' par(new=TRUE)
#' 
#' plot(density(pv$PV1), col='green', axes=FALSE, xlab=NA, ylab=NA, main=NA)
#' axis(side=4)
#' mtext(side = 4, line = 2.5, 'population density (green)')
#' 
#' par(op)
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
information = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='information')
}

#' @rdname information
expected_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='expected')
}

#' @rdname information
r_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='sim')
}

#' @rdname information
p_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='pmf')
}


theta_function = function(parms, items=NULL, booklet=NULL, which.draw=NULL, 
                          what=c('information','expected','sim','pmf'))
{
  what = match.arg(what)
  if(is.factor(items)) items = as.character(items)
  check_character(items,nullable=TRUE)
  check_string(booklet,name='booklet_id',nullable=TRUE)
  check_num(which.draw,nullable=TRUE)
  
  # data preparation
  # create fl(item_id,first,last), a, b
  
  if(inherits(parms,'data.frame'))
  {
    if(!is.null(items))
      parms = filter(parms,.data$item_id %in% items)
    
    if(!is.null(booklet))
      parms = filter(parms,.data$booklet_id == booklet)
    
    out = transform.df.parms(parms,'b')
    a = out$item_score
    b = out$b
    
    fl = out |>
      mutate(rn=row_number()) |>
      group_by(.data$item_id) |>
      summarise(first=as.integer(min(.data$rn)), last=as.integer(max(.data$rn))) |>
      ungroup() 
    

  } else if(inherits(parms,'prms'))
  {
    a = parms$inputs$ssIS$item_score
    b = parms$est$b
    if(is.matrix(b))
    {
      if(is.null(which.draw))
      {
        b = colMeans(b)
      } else
      {
        if(which.draw<1 || which.draw > nrow(b))
          stop('argument `which.draw` out of range')
        b = as.vector(b[which.draw,])
      }
    }
    
    fl = parms$inputs$ssI
    fl$item_id = as.character(fl$item_id)
    
    if(!is.null(items))
    {
      items = unique(items)
      suppressWarnings({fl = semi_join(fl, tibble(item_id=items),by='item_id')})
      if(nrow(fl) != length(items))
      {
        message('unknown items:')
        print(sort(setdiff(items,fl$item_id)))
        stop('Some items were not found, see output')
      }
    } else if(!is.null(booklet))
    {
      booklet = unique(booklet)
      design = parms$inputs$design
      if(length(intersect(booklet,design$booklet_id))<length(booklet))
      {
        stop('unknown booklet')
      }
      
      fl = design |>
        filter(.data$booklet_id %in% booklet) |>
        distinct(.data$item_id, .keep_all=TRUE)
    }  
    fl = arrange(fl,.data$first)
  } else stop_('parms must be a data.frame or an object of type `prms`') 
  rm(parms)  
  #output
  
  if(what=='information')
  {
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta)))
        stop('theta may not contain nan/NA values')  
      
      res = rep(0,length(theta))
      res[is.finite(theta)] = IJ_(b,a,fl$first, fl$last, theta[is.finite(theta)])$I
      # extremely large values overflow to NaN, recover as 0
      res[is.nan(res)] = 0
      res
    }
    class(out) = append('inf_func',class(out))
    
  } else if(what == 'expected')
  {
    max_score = sum(a[fl$last])
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta))) 
        stop('theta may not contain nan/NA values') 
      
      res = rep(0,length(theta))
      
      res[is.finite(theta)] = E_score(theta[is.finite(theta)],  
                                      b=b, a=a, 
                                      first=fl$first, last=fl$last)
      
      res[is.infinite(theta) & theta > 0] = max_score
      # extremely large values of theta overflow to NaN (small values undeflow to zero, which is fine)
      res[is.nan(res)] = max_score
      res
    }
    class(out) = append('exp_func',class(out))
  } else if(what=='sim')
  {
    out = function(theta)
    {
      res = rscore_item(theta,b=b,a=a,first = fl$first, last = fl$last)
      colnames(res) = fl$item_id
      res
    }
    class(out) = append('sim_func',class(out))
  } else if(what=='pmf')
  {
    out = function(theta)
    {
      res = pscore(theta,b=b,a=a,first = fl$first, last = fl$last)
      t(res)
    }
    class(out) = append('pmf_func',class(out))
  }
  
  out
}

print.inf_func = function(x,...) cat('Information function: I(theta)\n')
print.exp_func = function(x,...) cat('Conditional expected score function: E(X_i|theta)\n')
print.sim_func = function(x,...) cat('function to simulate item scores: (x_i1, ..., x_ip) ~ ENORM(theta)\n')
print.pmf_func = function(x,...) cat('Conditional score distribution function: P(x_+|theta)\n')


