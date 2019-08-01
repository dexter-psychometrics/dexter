
# to do: decide on one way to call merge_within_person

##########################################
#' Estimate abilities
#'
#' Computes estimates of ability for persons or booklets
#'
#' @param dataSrc Data source: a connection to a dexter database or a data.frame with columns or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by \code{\link{fit_enorm}} and containing
#' parameter estimates
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param design A data.frame with columns item_id and optionally booklet_id. If design is NULL
#' the score transformation table will be computed based on the test design 
#' that was used to calibrate the items. If the column booklet_id is not included, the score 
#' transformation table will be based on all items found in the design.
#' @param method   Maximum Likelihood (MLE) or Expected A posteriori (EAP) 
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
#' @param merge_within_person for persons who were administered multiple booklets, 
#' whether to provide just one abilty value (TRUE) or one per booklet(FALSE)
#' 
#' @return 
#' \describe{
#'   \item{ability}{a data.frame with columns: booklet_id, person_id, booklet_score, theta and optionally se (standard error) }
#'   \item{ability_tables}{a data.frame with columns: booklet_id, booklet_score, theta and optionally se (standard error)}
#' }
#' 
#' @details MLE estimates of ability will produce an NA for
#' the minimum (=0) or the maximum score on a booklet. If this is undesirable, 
#' we advise to use EAP with Jeffreys prior.
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
#' 
ability = function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP"), prior=c("normal", "Jeffreys"), 
                   use_draw=NULL, npv=500, mu=0, sigma=4, standard_errors=FALSE, merge_within_person=FALSE)
{
# to do: see if we can get factor level warnings with selected user data 
  check_dataSrc(dataSrc)
  check_parms(parms)
  
  method = match.arg(method)
  prior = match.arg(prior) 
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=TRUE, env=env, 
                           parms_check=parms$inputs$ssIS[,c('item_id','item_score')], merge_within_person = merge_within_person)
  
  if(nrow(respData$x) == 0) 
    stop('no response data to analyse')

  abl = ability_tables_(parms=parms, design = respData$design, method = method, prior=prior, use_draw = use_draw, 
                       npv=npv, mu=mu, sigma=sigma, standard_errors=standard_errors)
  
  respData$x %>% 
    inner_join(abl, by = c("booklet_id", "booklet_score")) %>% 
    select(suppressWarnings(one_of('booklet_id', 'person_id', 'booklet_score', 'theta', 'se'))) %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
}



#' @rdname ability
ability_tables = function(parms, design = NULL, method = c("MLE","EAP"), prior=c("normal", "Jeffreys"), 
                          use_draw = NULL, npv=500, mu=0, sigma=4, standard_errors = TRUE)
{
  method = match.arg(method)
  prior = match.arg(prior) 
  
  if(!is.null(design))
  {
    colnames(design) = tolower(colnames(design))
    
    if(! 'booklet_id' %in% colnames(design)) 
      design$booklet_id = 'all_items'
    
    if(! 'item_id' %in% colnames(design)) 
      stop('design must at least contain the column item_id')
    
    design = design %>%
      distinct(.data$booklet_id, .data$item_id) %>%
      mutate_if(is.factor, as.character)
  }
  
  if(method=='EAP' && prior=="normal")
  {
    check_num(npv, 'integer', .length=1, .min=1)
    check_num(mu, .length=1)
    check_num(sigma, .length=1)
  }
  check_num(use_draw, 'integer', .length=1, nullable=TRUE)
  
  if (method=="EAP")
  {
    if (sigma<0)
    {
      warning("Prior sd cannot be negative. Set to 4.")
      sigma = 4
    }
    if (npv<1)
    {
      warning("Number of plausible values must be positive. Set to 500")
      npv = 500
    }
  }
  
  
  ability_tables_(parms, design=design, method = method, 
                prior=prior, use_draw = use_draw, 
                npv=npv, mu=mu, sigma=sigma, standard_errors = standard_errors) %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
}


ability_tables_ = function(parms, design=NULL, method = c("MLE","EAP"), prior=c("normal", "Jeffreys"), 
                           use_draw = NULL, npv=500, mu=0, sigma=4, standard_errors = TRUE)
{

  method = match.arg(method)
  prior = match.arg(prior) 
  
  if (method=="EAP" && prior=="Jeffreys") 
    method="jEAP"
  
  if(is.null(design))
  {
    if(is.null(parms$inputs$design))
    {
      # old dexterMST version compatibility
      # dexterMST parms nor data has factors
      design = lapply(parms$inputs$bkList, function(bk) tibble(booklet_id=bk$booklet,item_id=bk$items)) %>% 
        bind_rows() %>% 
        inner_join(parms$inputs$ssI, by='item_id')  %>% 
        arrange(.data$booklet_id, .data$first)
    } else
    {
      design = parms$inputs$design
    }
  } else
  {
    # can have factor level warnings which we do not care about
    suppressWarnings({
      design = design %>%
        left_join(parms$inputs$ssI, by='item_id') %>% 
        arrange(.data$booklet_id, .data$first)
    })
    if(any(is.na(design$first))) 
      stop('some of the items in design are not present in the parms object')
    
  }

  a = parms$inputs$ssIS$item_score

  if(parms$inputs$method=="CML"){
    b = parms$est$b
  } else 
  {
    if(is.null(use_draw)) 
    {
      b = colMeans(parms$est$b)  
    } else 
    {
      b = parms$est$b[min(nrow(parms$est$b),use_draw),]	
    }   
  } 
  estimate = switch(method, 
     'MLE'  = function(.){ theta_MLE(b, a, .$first, .$last, se=standard_errors) }, 
     'EAP'  = function(.){ theta_EAP(b, a, .$first, .$last, npv=npv, mu=mu, sigma=sigma, se=standard_errors) }, 
     'jEAP' = function(.){ theta_jEAP(b, a, .$first, .$last, se=standard_errors) })
  
  # under the assumption that we always get theta's for the vector 0:max_test_score 
  design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = estimate(.)
      out = tibble(booklet_score=0:(length(est$theta)-1), theta = est$theta)
      if(standard_errors)
        out$se = est$se
      out
    }) %>%
    ungroup() 

}

