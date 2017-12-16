

##########################################
#' Generate plausible testscores
#'
#' Generate plausible i.e., posterior predictive sumscores on a set of items. 
#' A typical use of this function is to generate plausible scores on
#' a complete item bank when data is collected using an incomplete design
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates. If parms is given the function provides plausible scores conditional on the 
#' item parameters. These are considered known. If \code{parms} is \code{NULL}, Bayesian parameters are calculated from the datasrc
#' @param items vector of item_id's, this specifies the itemset to generate the testscores for. If \code{items} is \code{NULL} 
#' all items occurring in \code{dataSrc} are used.
#' @param covariates name or a vector of names of the variables to group the population, used to update the prior.
#' A covariate must be a discrete person covariate (e.g. not a float) that indicates nominal categories, e.g. gender or school
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param keep.observed In some cases, responses to one or more of the items have been observed. 
#' The user can choose to keep these observations or generate new ones.
#' @param nPS Number of plausible testscore to generate per person.
#' @return A data.frame with columns booklet_id, person_id, sumScore and nPS plausible scores
#' named PS1...PSn.
#'  
plausible_scores = function(dataSrc, predicate=NULL, parms=NULL, items=NULL, covariates=NULL, keep.observed=TRUE, nPS=1)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised = FALSE, extra_columns = covariates, env = env)
  
  # clean up items
  if(is.null(items)) items = respData$design$item_id
  if(inherits(items, 'data.frame')) items = items$item_id
  items = unique(items)
  
  # if no parms, we calculate parms from the full data source
  # and we opt for a Bayesian approach
  if(is.null(parms))
  {
    use_b_matrix = TRUE
    parms = fit_enorm_(respData, method='Bayes', nIterations=(nPS*5-1))
    b = parms$est$b[seq(4,5*nPS,by=5),]
    if (is.vector(b)) dim(b)=c(1,length(b)) 
  } else
  {
    if(length(setdiff(items, parms$inputs$ssI$item_id)) > 0)
    {
      message('The parms object does not contain parameter estimates for the following items')
      print(setdiff(items, parms$inputs$ssI$item_id))
      stop('Some items are without parameters')
    }
    
    if (parms$input$method=='CML')
    {
      use_b_matrix = FALSE
      b = parms$est$b
    }else
    {
      if ( nrow(parms$est$b)>=(nPS*5-1) )#are there enough rows?
      {
        use_b_matrix=TRUE
        b = parms$est$b[seq(4,5*nPS,by=5),]
        if (is.vector(b)) dim(b)=c(1,length(b)) 
      }else
      {
        use_b_matrix=FALSE
        b = colMeans(parms$est$b)
      }
    }
    # remove items from respdata for which we do not have parameters
    # otherwise they will cause an error in plausible value computation

    respData = semi_join(respData, parms$inputs$ssI, by='item_id', .recompute_sumscores = TRUE)
  }
  
  a = parms$inputs$ssIS$item_score
  
  # now we make plausible values using all responses we stil have
  pv = plausible_values_(respData, parms = parms, covariates = covariates, nPV = nPS, use_b_matrix = use_b_matrix)

  #save the design since respData may be mutilated below
  design = respData$design
  
  # remove responses that are not part of the specified itemset and recompute
  if(length(intersect(respData$design$item_id, items)) == 0)
  {
    pv = mutate(pv, sumScore = 0)
  } else
  {
    # since the sumscore in the plausible values is based on more responses than occur in the item set
    # we remove it in favor of the reduced testscores in the data
    # if persons made 0 items in the selected set they get a testscore of 0
    
    respData = semi_join(respData, tibble(item_id=items), by='item_id', .recompute_sumscores = TRUE)
    # summarise to one row per person
    respData = get_resp_data(respData, summarised = TRUE)
    pv = pv %>% 
      select(-.data$sumScore) %>%
      left_join(respData$x,  by=c("person_id", "booklet_id")) %>%
      mutate(sumScore = coalesce(.data$sumScore, 0L))
  }

  
  items = tibble(item_id = items) %>%
    inner_join(parms$inputs$ssI, by='item_id') %>%
    select(.data$item_id, .data$first, .data$last)
  
  if(keep.observed)
  {
    # we need to have a list of booklets containing first and last
    # for those items that were NOT in a booklet
    # this can not be based on the parms, since parms can have different booklets than data(and thus pv)
    # so we base it on design
    bkList = lapply(split(design, design$booklet_id), function(bk_items){ items %>% anti_join(bk_items, by='item_id') %>% arrange(.data$first)})


    pv = pv %>%
      group_by(.data$booklet_id) %>%
      do({
        bk = as.data.frame(bkList[[.$booklet_id[1]]])
        if(nrow(bk)==0)
        {
          mutate_at(.,vars(starts_with('PV')),`<-`, .$sumScore )
        } else
        {
          # we generate a score on the unobserved items and add to the score on the observed items
          cnt = counter$new()
          mutate_at(.,vars(starts_with('PV')), rscore, b=b, a=a, first=bk$first, last=bk$last, cntr=cnt, use_b_matrix = use_b_matrix) %>%
            mutate_at(vars(starts_with('PV')), `+`, .$sumScore)
        }
      }) %>%
      ungroup()
  } else
  {
    cnt = counter$new()
    pv = pv %>%
      mutate_at(vars(starts_with('PV')), rscore, b=b, a=a, first=items$first,last=items$last,cntr=cnt, use_b_matrix = use_b_matrix)
  }
  
  pv %>%
    rename_at(vars(starts_with('PV')), function(x) sub('PV','PS', x)) %>%
    select(-.data$sumScore)
}


