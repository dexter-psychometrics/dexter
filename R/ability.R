
##########################################
#' Estimate abilities
#'
#' Computes estimates of ability for persons or booklets
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
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
#' @param standard_errors If true standard-errors are produced.
#' @param asOPLM Report abilities on a scale defined by the OPLM normalization. Only when no values in parms have been fixed
#' 
#' @return 
#' \describe{
#'   \item{ability}{a data.frame with columns: booklet_id, person_id, sumScore, theta and optionally se (standard error) }
#'   \item{ability_tables}{a data.frame with columns: booklet_id, sumScore, theta and optionally se (standard error)}
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
#' plot(bb$sumScore, bb$theta, xlab="test-score", ylab="ability est.", pch=19, cex=0.7)
#' points(aa$sumScore, aa$theta, col="red", pch=19, cex=0.7)
#' points(aa$sumScore, cc$theta, col="green", pch=19, cex=0.7)
#' legend("topleft", legend = c("EAP normal prior", "EAP Jeffreys prior", "MLE"), bty = "n",
#'         lwd = 1, cex = 0.7, col = c("black", "green", "red"), lty=c(0,0,0), pch = c(19,19,19))
#' 
#' close_project(db)
#' }
#' 
#' 
ability = function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP"), prior=c("normal", "Jeffreys"), use_draw=NULL, 
                    npv=500, mu=0, sigma=4, standard_errors=FALSE,  asOPLM=TRUE){

  check_arg(dataSrc, 'dataSrc')
  check_arg(parms, 'prms')
  
  method <- match.arg(method)
  prior = match.arg(prior) 
  qtpredicate = eval(substitute(quote(predicate)))
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=caller_env())
  
  if(nrow(respData$x)==0) stop('no data to analyse')
  
  # check if all items and scores are known
  if(nrow(anti_join(respData$design, parms$inputs$ssI, by=c('item_id'))) > 0)
  {
    stop('some of your items are without parameters')
  } 
  if(nrow(anti_join(respData$x, parms$inputs$ssIS, by=c('item_id','item_score'))) > 0)
  {
    stop('Some item_scores in your data are not present in your parameters.')
  } 
  # now we can summarise
  respData = get_resp_data(respData, summarised = TRUE)

  abl = ability_tables(parms=parms, design = respData$design, method = method, prior=prior, use_draw = use_draw, 
                       npv=npv, mu=mu, sigma=sigma, standard_errors=standard_errors, asOPLM=asOPLM)
  
  return(respData$x %>% 
           inner_join(abl, by = c("booklet_id", "sumScore")) %>% 
           select(suppressWarnings(one_of('booklet_id', 'person_id', 'sumScore', 'theta', 'se'))) %>%
           as.data.frame())
}



#' @rdname ability
ability_tables = function(parms, design = NULL, method = c("MLE","EAP"), prior=c("normal", "Jeffreys"), use_draw = NULL, 
                          npv=500, mu=0, sigma=4, standard_errors = TRUE, asOPLM=TRUE) #smooth=FALSE,
{

  method = match.arg(method)
  prior = match.arg(prior) 
  if (method=="EAP" && prior=="Jeffreys") method="jEAP"
  if(method=='EAP' && prior=="normal")
  {
    check_arg(npv, 'integer', .length=1)
    check_arg(mu, 'numeric', .length=1)
    check_arg(sigma, 'numeric', .length=1)
  }
  check_arg(asOPLM, 'logical', .length=1)
  check_arg(standard_errors, 'logical', .length=1)
  check_arg(use_draw, 'integer', .length=1, nullable=TRUE)
  
  if (asOPLM && !parms$inputs$has_fixed_parms)
  {
    ff = toOPLM(parms$inputs$ssIS$item_score, parms$est$b, parms$inputs$ssI$first, parms$inputs$ssI$last)
    if (parms$inputs$method=="CML"){
      parms$est$b = toDexter(parms$est$beta.cml, ff$a, ff$first, ff$last, re_normalize = FALSE)$est$b
    }else
    {
      for (i in 1:nrow(parms$est$b))
      {
        parms$est$b[i,] = toDexter(parms$est$beta.cml[i,], ff$a, ff$first, ff$last, re_normalize = FALSE)$est$b
      }
    }
  }
  
  ## Check validity of input
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
  
  if(is.null(design))
  {
    if(is.null(parms$inputs$design))
    {
      # old dexterMST version compatibility
      design = lapply(parms$inputs$bkList, function(bk) tibble(booklet_id=bk$booklet,item_id=bk$items)) %>% bind_rows()
    } else
    {
      design = select(parms$inputs$design, .data$booklet_id, .data$item_id)
    }
  } else 
  {
    # clean up the design if necessary
    colnames(design) = tolower(colnames(design))
    if(! 'booklet_id' %in% colnames(design)) design$booklet_id = 'all_items'
    if(! 'item_id' %in% colnames(design)) stop('design must at least contain the column item_id')
    design = design %>%
      select(.data$booklet_id, .data$item_id) %>%
      mutate_if(is.factor, as.character)
  }
  # should we make (booklet_id, item_id) unique or is the user allowed footshooting here?
  
  design = design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  if(parms$inputs$method=="CML"){
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
  } else {
    a = parms$est$a
    if(is.null(use_draw)) {
      b = colMeans(parms$est$b)  
    } else {
      b = parms$est$b[min(nrow(parms$est$b),use_draw),]	
    }   
  } 
  # add b to ssIS
  parms$inputs$ssIS$b = b
  estimate = switch(method, 
                    'MLE'  = function(.){ theta_MLE(b, a, .$first, .$last, se=standard_errors) }, 
                    'EAP'  = function(.){ theta_EAP(b, a, .$first, .$last, npv=npv, mu=mu, sigma=sigma, se=standard_errors) }, 
                    'jEAP' = function(.){ theta_jEAP(b, a, .$first, .$last, se=standard_errors) })
  
  # under the assumption that we always get theta's for the vector 0:max_test_score 
  design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = estimate(.)
      
      if(standard_errors)
      {
        tibble(booklet_id = rep_len(.$booklet_id, length(est$theta)),
               sumScore=0:(length(est$theta)-1),
               theta = est$theta,
               se = est$se)
      } else
      {
        tibble(booklet_id = rep_len(.$booklet_id, length(est$theta)),
                   sumScore=0:(length(est$theta)-1),
                   theta = est$theta)
      }
    }) %>%
    ungroup() %>%
    as.data.frame()
}



## Experimental: produces for unweighted scores
theta_tables = function(parms, design = NULL, method = c("MLE","EAP"), use_draw = NULL, 
                        npv=500, mu=0, sigma=4, asOPLM=TRUE) 
{
  method = match.arg(method)
  if ((asOPLM)&(!parms$inputs$has_fixed_parms))
  {
    ff = toOPLM(parms$inputs$ssIS$item_score, parms$est$b, parms$inputs$ssI$first, parms$inputs$ssI$last)
    if (parms$inputs$method=="CML"){
      parms$est$b = toDexter(parms$est$beta.cml, ff$a, ff$first, ff$last, re_normalize = FALSE)$est$b
    }else
    {
      for (i in 1:nrow(parms$est$b))
      {
        parms$est$b[i,] = toDexter(parms$est$beta.cml[i,], ff$a, ff$first, ff$last, re_normalize = FALSE)$est$b
      }
    }
  }
  
  if(is.null(design))
  {
    design = lapply(parms$inputs$bkList, function(bk) tibble(booklet_id=bk$booklet,item_id=bk$items)) %>% bind_rows()
  } else 
  {
    # clean up the design if necessary
    colnames(design) = tolower(colnames(design))
    if(! 'booklet_id' %in% colnames(design)) design$booklet_id = 'all_items'
    if(! 'item_id' %in% colnames(design)) stop('design must at least contain the column item_id')
    design = design %>%
      select(.data$booklet_id, .data$item_id) %>%
      mutate_if(is.factor, as.character)
  }
  # should we make (booklet_id, item_id) unique or is the user allowed footshooting here?
  
  design = design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  if(parms$inputs$method=="CML"){
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
  } else {
    a = parms$est$a
    if(is.null(use_draw)) {
      b = colMeans(parms$est$b)  
    } else {
      b = parms$est$b[min(nrow(parms$est$b),use_draw),]	
    }   
  } 
  # add b to ssIS
  parms$inputs$ssIS$b = b
  A = unlist(lapply(parms$inputs$ssI$nCat, function(x)c(0:(x-1))),use.names = F) ## a for unweighted scores
  estimate = switch(method, 
                    'MLE'  = function(.){ theta_aA(b, a, A, .$first, .$last) }, 
                    'EAP'  = function(.){ theta_EAP(b, a, .$first, .$last, npv=npv, mu=mu, sigma=sigma, se=TRUE, A) })
 
  # under the assumption that we always get theta's for the vector 0:max_test_score 
  design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = estimate(.)
      if(method=="EAP")
      {
        tibble(booklet_id = rep_len(.$booklet_id, length(est$theta)),
               sumScore=0:(length(est$theta)-1),
               theta = est$theta,
               se = est$se)
      } else
      {
        tibble(booklet_id = rep_len(.$booklet_id, length(est$theta)),
               sumScore=0:(length(est$theta)-1),
               theta = est$theta)
      }
    }) %>%
    ungroup()
}


