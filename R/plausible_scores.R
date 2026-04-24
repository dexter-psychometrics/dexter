
#' Draw plausible test scores
#'
#' Draw plausible, i.e. posterior predictive sumscores on a set of items. 
#' 
#' A typical use of this function is to generate plausible scores on
#' a complete item bank when data is collected using an incomplete design
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates. If parms is given the function provides plausible scores conditional on the 
#' item parameters. These are considered known. If \code{parms} is \code{NULL}, Bayesian parameters are calculated from the dataSrc
#' @param parms_draw when the item parameters are estimated Bayesianly (see: \code{\link{fit_enorm}}), 
#' parms_draw specifies whether to use a sample(a different item parameter draw for each plausible values draw) or the posterior mean
#' of the item draws. Alternatively, it can be an integer specifying a specific draw. Ignored when parms is not estimated Bayesianly.
#' @param link_error only for cml, whether to use the maximum likelihood estimates or random draws based on the covariance matrix, 
#' thus including estimation/linking error in the pv draws.
#' @param items vector of item_id's, this specifies the itemset to generate the testscores for. If \code{items} is \code{NULL} 
#' all items occurring in \code{dataSrc} are used.
#' @param covariates name or a vector of names of the variables to group the population, used to update the prior.
#' A covariate must be a discrete person covariate that indicates nominal categories, e.g. gender or school
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param keep.observed If responses to one or more of the items have been observed,
#' the user can choose to keep these observations or generate new ones. 
#' @param nPS Number of plausible testscores to generate per person.
#' @param prior_dist use a normal prior for the plausible values or a mixture of two normals. 
#' A mixture is only possible when there are no covariates.
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible scores are generated per person (TRUE) or per booklet (FALSE)
#' @param by_item return scores per item instead of sumscores
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPS plausible scores
#' named PS1...PSn.
#'  
plausible_scores = function(dataSrc, parms=NULL, predicate=NULL, items=NULL, parms_draw = c('sample','average'),
                            link_error=FALSE,
                            covariates=NULL, nPS=1, prior_dist = c("normal", "mixture"),
                            keep.observed=TRUE, by_item=FALSE, merge_within_persons=FALSE)  
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  check_dataSrc(dataSrc)
  parms_draw = match.arg(parms_draw)
  prior_dist = match.arg(prior_dist)
  
  df_info = get_datatype_info(dataSrc, columns = c('booklet_id','item_id','person_id',covariates))
  
  pb = get_prog_bar(nsteps=if(is.null(parms)) 130 else 100, 
                    retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env,
                           merge_within_persons=merge_within_persons)
  
  if(is.null(items))
  {
    items = levels(respData$design$item_id)
  } else if(inherits(items,'data.frame'))
  {
    items = as.character(unique(items$item_id))
  } else
  {
    items = as.character(unique(items))
  }
  
  # if there are no params, all of items must be in data
  # if there are params, all of items must be in params
  
  if(is.null(parms) && !all(items %in% levels(respData$design$item_id)))
  {
    stop_("`items` contains item_id's not found in the data, you must either provide parameters reparately or ",
          "specify only items present in your data")
  } else if(!is.null(parms))
  {
    if(inherits(parms,'data.frame')) parms_items = as.character(unique(parms$item_id))
    else parms_items = unique(coef(parms)$item_id)
    
    if(!all(items %in% parms_items))
      stop_("`items` contains item_id's not found in the parameters")
  }
  
  # generate plausible values and params
  res = plausible_values_(respData, parms=parms, covariates=covariates, 
                          nPV=nPS, parms_draw = parms_draw, link_error=link_error,
                          prior_dist = prior_dist)
  
  parms = res$parms
  pv = res$pv

  fl = parms$items |>
    filter(.data$item_id %in% items) |>
    select('item_id','first0','last0')
  
  if(keep.observed)
  {
    # unfortunately cannot get around a sort
    xist = respData$x |>
      inner_join(select(fl,'item_id','first0'), by='item_id') |>
      inner_join(tibble(person_id = pv$person_id, booklet_id=pv$booklet_id, pv_indx=0:(nrow(pv)-1L)), 
                 by=c('booklet_id', 'person_id')) |>
      select(person_id = 'pv_indx',item_first = 'first0','item_score') |>
      arrange(.data$person_id, .data$item_first)
    
  } else
  {
    xist = tibble(person_id=-1L,item_first=-1L,item_score=-1L)
  }

  ps = sample_scores(theta= as.matrix(select(pv,matches('^PV\\d+$'))), parms$b, parms$a, fl$first0, fl$last0, 
                           by_item = by_item, item_long=TRUE,
                           existing_scores = xist)
  
  colnames(ps) = paste0('PS',1:ncol(ps))
  
  if(by_item)
  {
    indx = rep(1:nrow(pv), each=nrow(fl))
    x = select(pv, -'booklet_score', -matches('^PV\\d+$'))[indx,]
    x$item_id = rep(fl$item_id,nrow(pv))

  } else
  {
    x = select(pv, -'booklet_score', -matches('^PV\\d+$')) 
  }

  bind_cols(x, ps) |>
      df_format(df_info)
}
