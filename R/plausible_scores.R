
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
                          nPV=nPS, parms_draw = parms_draw, 
                          prior_dist = prior_dist)
  
  parms = res$parms
  pv = res$pv
  
  items = sort(factor(items,levels=levels(parms$items$item_id)))
  
  # 0-indexed for use in c functions
  fl = parms$items |>
    filter(.data$item_id %in% items) |>
    select('item_id','first0','last0')
  
  multiple_b = !is.null(dim(parms$b)) && ncol(parms$b)>1
  
  a = parms$a
  if(multiple_b)
  {
    b = t(parms$b)
    bstep = as.integer((ncol(b)-1)/max(nPS-1,1))
  } else
  {
    b = matrix(parms$b,ncol=1)
    bstep = 0L
  }
  
  if(by_item)
  {
    if(merge_within_persons)
    {
      pp = tibble(old_person_id=pv$person_id,booklet_id=pv$booklet_id, person_id=1:nrow(pv))
      pv$person_id = pp$person_id
      
      respData$x = respData$x |> 
        rename(old_person_id='person_id') |>
        inner_join(pp,by=c('booklet_id','old_person_id')) |>
        select(-'old_person_id')
    }
    
    pv = droplevels(pv)
    if(is.unsorted(pv$person_id)) pv = arrange(pv,.data$person_id)

    nit = length(items)
    np = nrow(pv)
    
    x = tibble(person_id=rep(pv$person_id,nit),item_id=rep(items,each=np)) 
      
    b_index = 1L
    for(i in 1:nPS)
    {
        
      x[[sprintf('PS%i',i)]] = as.integer(sampleNRM_itemC(pv[[sprintf('PV%i',i)]], b[,b_index], a, fl$first0, fl$last0))
      b_index = b_index + bstep
    }
    
    if(keep.observed && any(respData$design$item_id %in% items))
    {
      m = get_resp_matrix(filter(respData$x, .data$item_id %in% items ))
      if(ncol(m) == length(items))
      {
        w = which(!is.na(m))
        m = as.integer(m[w])
      }
      else
      {
        # not all items are in responses
        w = apply(m,2,\(mcol) unname(which(!is.na(mcol))),simplify=FALSE)
        w = mapply(w,0:(length(w)-1), FUN= \(a,b) a+b*np, SIMPLIFY=FALSE )

        ii = left_join(tibble(item_id=items), tibble(item_id=colnames(m), present=TRUE), by='item_id') |>
          mutate(present=coalesce(.data$present,FALSE)) |>
          arrange(.data$item_id) |>
          mutate(cnt=cumsum(as.integer(!.data$present))*np) |>
          filter(.data$present)
          
        m = as.integer(m[unlist(w)])
          
        w = unlist(mapply(ii$cnt, w, FUN='+',SIMPLIFY=FALSE))
      }
      
      for(i in 1:nPS)
      {
        x[[sprintf('PS%i',i)]][w] = m
      }
    }
    
    if(!is.null(covariates))
      pv = inner_join(select(pv,all_of(c('person_id',covariates))), x, by='person_id')
    else
      pv = x
    
    if(merge_within_persons)
    {
      pv = inner_join(pv,pp,by='person_id') |>
        select(-'person_id') |>
        mutate(person_id='old_person_id') 
    }
    
  } else
  {
    if(keep.observed && any(respData$design$item_id %in% items))
    {
      # keep track of sumscore on selected items
      
      respData = semi_join_rd(respData, tibble(item_id=items), by='item_id', .recompute_sumscores = TRUE)
      respData = get_resp_data(respData, summarised = TRUE, protect_x = !is_db(dataSrc))
      
      pv = pv |> 
        select(-'booklet_score') |>
        left_join(respData$x,  by=c("person_id", "booklet_id")) |>
        mutate(booklet_score = coalesce(.data$booklet_score, 0L))
      
      pv = lapply(split(pv, pv$booklet_id), function(pvbk)
      {
        bk = pvbk$booklet_id[1]
        
        fl_bk = fl |>
          anti_join(filter(respData$design, .data$booklet_id == bk), by='item_id')
        
        #nothing to augment case
        if(nrow(fl_bk) == 0)
        {
          for(pn in sprintf('PV%i',1:nPS)) pvbk[[pn]] = pvbk$booklet_score
        } else
        {
          b_index = 1L
          for(pn in sprintf('PV%i',1:nPS))
          {
            pvbk[[pn]] = sampleNRM_testC(pvbk[[pn]], b[,b_index], a, fl_bk$first0, fl_bk$last0)[,1,drop=TRUE] + pvbk$booklet_score
            b_index = b_index + bstep
          }    
        }
        pvbk
      }) |>
        bind_rows()
      
    } else
    {
      b_index = 1L
      
      for(pn in sprintf('PV%i',1:nPS))
      {
        pv[[pn]] = sampleNRM_testC(pv[[pn]], b[,b_index], a, fl$first0, fl$last0)[,1,drop=TRUE]
        b_index = b_index + bstep
      }
    }
  }

  pv |>
    select(-any_of('booklet_score')) |>
    rename_with(gsub, pattern='^PV(?=\\d+$)',replacement='PS', perl=TRUE)  |>
    df_format(df_info)
}
