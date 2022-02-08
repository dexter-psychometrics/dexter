


##########################################
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
#' @param use_draw When parms is Bayesian, use_draw is 
#'                 the index of the posterior sample of the item 
#'                 parameters that will be used for generating plausible values. 
#'                 If use_draw=NULL, a posterior mean is used. 
#'                 If outside range, the last iteration will be used. 
#' @param npv Number of plausible values sampled to calculate EAP with normal prior
#' @param mu Mean of the normal prior
#' @param sigma Standard deviation of the normal prior
#' @param standard_errors If true standard-errors are produced
#' @param merge_within_persons for persons who were administered multiple booklets, 
#' whether to provide just one ability value (TRUE) or one per booklet(FALSE)
#' 
#' @return 
#' \describe{
#'   \item{ability}{a data.frame with columns: booklet_id, person_id, booklet_score, theta and optionally se (standard error) }
#'   \item{ability_tables}{a data.frame with columns: booklet_id, booklet_score, theta and optionally se (standard error)}
#' }
#' 
#' @details MLE estimates of ability will produce an NA for
#' the minimum (=0) or the maximum score on a booklet. If this is undesirable, 
#' we advise to use WLE. The WLE was proposed by Warm (1989) to reduce bias in the MLE and is also known
#' as the Warm estimator.
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' f = fit_enorm(db)
#' aa = ability_tables(f,method="MLE",standard_errors=FALSE)
#' bb = ability_tables(f,method="EAP",standard_errors=FALSE)
#' cc = ability_tables(f,method="EAP",prior="Jeffreys", standard_errors=FALSE)
#' plot(bb$booklet_score, bb$theta, xlab="test-score", ylab="ability est.", pch=19, cex=0.7)
#' points(aa$booklet_score, aa$theta, col="red", pch=19, cex=0.7)
#' points(aa$booklet_score, cc$theta, col="green", pch=19, cex=0.7)
#' legend("topleft", legend = c("EAP normal prior", "EAP Jeffreys prior", "MLE"), bty = "n",
#'         lwd = 1, cex = 0.7, col = c("black", "green", "red"), lty=c(0,0,0), pch = c(19,19,19))
#' 
#' close_project(db)
#' }
#' 
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory. 
#' Psychometrika, 54(3), 427-450. 
#' 
ability = function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP","WLE"), prior=c("normal", "Jeffreys"), 
                   use_draw=NULL, mu=0, sigma=4, standard_errors=FALSE, merge_within_persons=FALSE)
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
  

  abl = ability_tables(parms=parms, design = respData$design, method = method, prior=prior, use_draw = use_draw, 
                       mu=mu, sigma=sigma, standard_errors=standard_errors)
  abl$booklet_id = ffactor(abl$booklet_id, levels = levels(respData$design$booklet_id))
  respData$x %>% 
    inner_join(abl, by = c("booklet_id", "booklet_score")) %>% 
    select(suppressWarnings(one_of('booklet_id', 'person_id', 'booklet_score', 'theta', 'se'))) %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}



#' @rdname ability
ability_tables = function(parms, design = NULL, method = c("MLE","EAP","WLE"), prior=c("normal", "Jeffreys"), 
                          use_draw = NULL, mu=0, sigma=4, standard_errors = TRUE)
{
  method = match.arg(method)
  prior = match.arg(prior) 
  
  if(method=='EAP' && prior=="normal")
  {
    check_num(mu, .length=1)
    check_num(sigma, .length=1, .min=0)
    check_num(use_draw, 'integer', .length=1, nullable=TRUE)
  }
  
  if (method=="EAP" && prior=="Jeffreys") method="jEAP"
  
  simple_parms = simplify_parms(parms, design, use_draw, collapse_b=TRUE) 
  b = simple_parms$b
  a = simple_parms$a
  
  estimate = switch(method, 
                    'MLE'  = function(.){ theta_MLE(b, a, .$first, .$last, se=standard_errors) }, 
                    #'EAP'  = function(.){ theta_EAP(b, a, .$first, .$last, npv=npv, mu=mu, sigma=sigma, se=standard_errors) }, 
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
