


#' Estimate abilities
#'
#' Computes estimates of ability for persons or for booklet scores
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and, 
#' beta
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
  
  if(inherits(parms,'enorm') || inherits(parms,'prms'))
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


