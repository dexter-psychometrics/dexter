


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
  
  respData$x %>% 
    inner_join(abl, by = c("booklet_id", "booklet_score")) %>% 
    select(any_of(c('booklet_id', 'person_id', 'booklet_score', 'theta', 'se'))) %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}



#' @rdname ability
ability_tables = function(parms, design = NULL, method = c("MLE","EAP","WLE"), prior=c("normal", "Jeffreys"), 
                          parms_draw = 'average', mu=0, sigma=4)
{
  standard_errors = TRUE
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
                    'MLE'  = function(.){ theta_MLE(b, a, .$first, .$last, se=standard_errors) }, 
                    'EAP'  = function(.){ theta_EAP_GH(b, a, .$first, .$last, mu=mu, sigma=sigma) },
                    'jEAP' = function(.){ theta_jEAP(b, a, .$first, .$last, se=standard_errors) },
                    'WLE' = function(.){ theta_WLE(b, a, .$first, .$last, se=standard_errors) })
  
  
  # under the assumption that we always get theta's for the vector 0:max_test_score 
  simple_parms$design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = estimate(.)
      out = tibble(booklet_score=0:(length(est$theta)-1), theta = est$theta)
      if(standard_errors)
        out$se = est$se
      out
    }) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}
