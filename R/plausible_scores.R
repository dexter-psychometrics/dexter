

##########################################
#' Generate plausible test scores
#'
#' Generate plausible i.e., posterior predictive sumscores on a set of items. 
#' A typical use of this function is to generate plausible scores on
#' a complete item bank when data is collected using an incomplete design
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
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
#' @param keep.observed If responses to one or more of the items have been observed,
#' the user can choose to keep these observations or generate new ones. 
#' @param nPS Number of plausible testscores to generate per person.
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible scores are generated per person (TRUE) or per booklet (FALSE)
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPS plausible scores
#' named PS1...PSn.
#'  
plausible_scores = function(dataSrc, parms=NULL, predicate=NULL, items=NULL, 
                            covariates=NULL, keep.observed=TRUE, nPS=1,merge_within_persons=FALSE)  
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  check_dataSrc(dataSrc)

  plausible_scores_(dataSrc, parms=parms, qtpredicate=qtpredicate, items=items, covariates=covariates, 
                    keep.observed=keep.observed, nPS=nPS, env=env,
                    merge_within_persons=merge_within_persons) %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}

# to do: this has become a messy function
plausible_scores_ = function(dataSrc, parms=NULL, qtpredicate=NULL, items=NULL, 
                             covariates=NULL, keep.observed=TRUE, nPS=1, env=NULL,
                             merge_within_persons=FALSE)
{
  if (is.null(env))
    env = caller_env()
  
  from = Gibbs.settings$from.ps 
  step = Gibbs.settings$step.ps # from and step are burnin and thinning
  keep.which = seq(from,(from-step)*(from>step)+step*nPS,by=step)
  nPS.needed = max(keep.which) # Given from and step, this many PS must be generated
  
  pb = get_prog_bar(nsteps=if(is.null(parms)) 120 else 100, 
                    retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  if(is.null(parms))
  {
    pcheck=NULL
  } else if(inherits(parms,'data.frame'))
  {
    parms = transform.df.parms(parms,'b', TRUE)
    pcheck = parms[,c('item_id','item_score')]
  } else
  {
    pcheck = parms$inputs$ssIS[,c('item_id','item_score')]
  }
  
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised = FALSE, 
                           extra_columns = covariates, env = env, 
                           parms_check=pcheck,
                           merge_within_persons=merge_within_persons)
  
  use_b_matrix = FALSE
  
  # if no parms, we calculate parms from the full data source
  # and we opt for a Bayesian approach
  if(is.null(parms))
  {
    nIter.enorm  = Gibbs.settings$from.pv + Gibbs.settings$step.pv*(nPS.needed-1)
    pb$new_area(20)
    parms = fit_enorm_(respData, method='Bayes', nDraws = nIter.enorm)
    
  } else if(inherits(parms,'prms') && parms$inputs$method != 'CML')
  {
    # if user parms was Bayesian we do some checks
    if (nrow(parms$est$b) < nPS.needed) # are there enough samples ?
    {
      stop(paste("Not enough posterior samples in fit_enorm for", 
              nPS, "plausible scores. Use at least", nPS.needed, "samples in fit_enorm"))

    } else
    {
      use_b_matrix = TRUE
    }
  }  
  
  # now we make plausible values using all responses we have
  pb$new_area(80) 
  pv = plausible_values_(respData, parms = parms, covariates = covariates, nPV = nPS.needed)
  pb$close_area()
  
  # clean up items
  if(is.null(items))
  {
    items = unique(respData$design$item_id)
  }
  else 
  {
    if(inherits(items, 'data.frame'))
      items = items$item_id
    
    items = unique(as.character(items))
  }
  
  simple_parms = simplify_parms(parms, collapse_b = !use_b_matrix, design = tibble(item_id=items))
  items = select(simple_parms$design, -.data$booklet_id) 
  a = simple_parms$a
  b = simple_parms$b
  items$item_id = re_factor_item_id(respData, items$item_id)
  levels(respData$x$item_id) = levels(respData$design$item_id) = levels(items$item_id)
  
  #save the full design since respData is be mutilated below
  design = respData$design
  
  # remove responses that are not part of the specified itemset and recompute
  if(length(intersect(respData$design$item_id, items$item_id)) == 0)
  {
    pv = mutate(pv, booklet_score = 0L)
  } else
  {
    # since the sumscore in the plausible values is (possibly) based on more responses than occur in the item set
    # we remove it in favor of the reduced testscores in the data
    # if persons made 0 items in the selected set they get a testscore of 0 through the left join below
    
    respData = semi_join(respData, items, by='item_id', .recompute_sumscores = TRUE)
    
    # summarise to one row per person
    respData = get_resp_data(respData, summarised = TRUE, protect_x=FALSE)
    
    pv = pv %>% 
      select(-.data$booklet_score) %>%
      left_join(respData$x,  by=c("person_id", "booklet_id")) %>%
      mutate(booklet_score = coalesce(.data$booklet_score, 0L))
  }
  pb$tick()
  if(keep.observed)
  {
    # we need to have a list of booklets containing first and last
    # for those items that were NOT in a booklet
    # this can not be based on the parms, since parms can have different booklets than data(and thus pv)
    # so we base it on design
    bkList = lapply(split(design, design$booklet_id), 
                    function(bk_items){ items %>% anti_join(bk_items, by='item_id') %>% arrange(.data$first)})
    
    # to do: message when nothing new is generated 
    pv = pv %>%
      group_by(.data$booklet_id) %>%
      do({
        bk = as.data.frame(bkList[[.$booklet_id[1]]])
        if(NROW(bk)==0)
        {
          mutate_at(.,vars(starts_with('PV')),`<-`, .$booklet_score )
        } else
        {
          # we generate a score on the unobserved items and add to the score on the observed items
          cntr = (function(){i=0L; function(){i<<-i+1L; i}})()
          mutate_at(.,vars(starts_with('PV')), rscore, b=b, a=a, first=bk$first, last=bk$last, cntr=cntr, use_b_matrix = use_b_matrix) %>%
            mutate_at(vars(starts_with('PV')), `+`, .$booklet_score)
        }
      }) %>%
      ungroup()
  } else
  {
    cnt =  (function(){i=0L;function(){i<<-i+1L;i}})()
    pv = pv %>%
      mutate_at(vars(starts_with('PV')), rscore, b=b, a=a, first=items$first,last=items$last,cntr=cnt, use_b_matrix = use_b_matrix)
  }
  
  pv = select(pv, all_of(covariates), .data$booklet_id,.data$person_id, matches('^PV\\d+$')) %>%
    rename_with(gsub,pattern='^PV(?=\\d+$)',replacement='PS', perl=TRUE) 

}



