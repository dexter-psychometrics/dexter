
# to do: example
##########################################
#' Profile analysis
#'
#' Expected and observed domain scores, conditional on the test score, per person or test score. 
#' Domains are specified as categories of items using item_properties.
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score, 
#' an arbitrarily named column containing an item property and optionally booklet_id
#' @param parms An object returned by \code{\link{fit_enorm}} 
#' @param predicate An optional expression to subset data in dataSrc, if NULL all data is used
#' @param item_property the name of the item property used to define the domains. If \code{dataSrc} is a dexter db then the
#' item_property must match a known item property. If datasrc is a data.frame, item_property must be equal to
#'  one of its column names. For profile_tables item_property must match a column name in \code{domains}.
#' @param design data.frame with columns item_id and optionally booklet_id
#' @param domains data.frame with column item_id and a column with name equal to \code{item_property} 
#' 
#' @return 
#' \describe{
#' \item{profiles}{a data.frame with columns person_id, booklet_id, booklet_score, 
#' <item_property>, domain_score, expected_domain_score}
#' \item{profile_tables}{a data.frame with columns booklet_id, booklet_score, 
#' <item_property>, expected_domain_score }
#' }
#' @details 
#' When using a unidimensional IRT Model like the extended nominal response model in 
#' dexter (see: \code{\link{fit_enorm}}), the model is as a rule to simple to catch all the relevant dimensions in a test.
#' Nevertheless, a simple model is quite useful in practice. Profile analysis can complement the model
#' in this case by indicating how a test-taker, conditional on her/his test score, 
#' performs on a number of pre-specified domains, e.g. in case of a mathematics test 
#' the domains could be numbers, algebra and geometry or in case of a digital test the domains could be animated versus
#'  non-animated items. This can be done by comparing the achieved score on a domain with the expected score, given the test score.
#' 
#' 
#' @references 
#' Verhelst, N. D. (2012). Profile analysis: a closer look at the PISA 2000 reading data. 
#' Scandinavian Journal of Educational Research, 56 (3), 315-332.
#' 
#' 
#
# to do: allow rim as parms
profiles = function(dataSrc, parms, item_property, predicate=NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  check_dataSrc(dataSrc)
  check_parms(parms)
  check_string(item_property)


  if(inherits(dataSrc, 'DBIconnection'))
    item_property = dbValid_colnames(item_property)
  
  # to do: is there an error if extra column is not available?
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env, 
                           extra_columns = item_property, 
                           parms_check=parms$inputs$ssIS[,c('item_id','item_score')])
  
  if(!is.integer(respData$x[[item_property]]))
    respData$x[[item_property]] = ffactor(respData$x[[item_property]])

  # check for valid item property
  if(!inherits(dataSrc, 'DBIConnection') || !item_property %in% dbListFields(dataSrc,'dxitems'))
  {
    idom = distinct(respData$x, .data$item_id, .data[[item_property]])
    
    if(nrow(idom) > nrow(distinct(idom, .data$item_id)))
      stop('Each item must have a unique item property')
  } 
  
  
  if(inherits(dataSrc,'DBIConnection'))
  {
    domains = dbGetQuery(dataSrc, paste("SELECT item_id,", item_property, "FROM dxItems;"))
    domains$item_id = ffactor(domains$item_id, levels = levels(respData$x$item_id))
    domains = filter(domains, !is.na(.data$item_id))
    if(!is.integer(domains[[item_property]]))
      domains[[item_property]] = ffactor(domains[[item_property]])
  } else
  {
    domains = distinct(respData$x, .data$item_id, .data[[item_property]])
    
    if(nrow(domains) > n_distinct(domains$item_id))
      stop("some items belong to multiple domains, this is not allowed")
  }
  
  design = respData$design
  
  respData = respData %>%
    polytomize(item_property, protect_x=!inherits(dataSrc,'DBIConnection')) 
  
  out = respData$x %>%
    inner_join(
        profile_tables_(parms = parms, design = design,
                       domains = domains,
                       item_property = item_property),
        by = c('booklet_id','booklet_score',item_id = item_property)) %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
  
  colnames(out)[colnames(out)=='item_id'] = item_property
  colnames(out)[colnames(out)=='item_score'] = 'domain_score'
  
  out
}
  
#' @rdname profiles
profile_tables = function(parms, domains, item_property, design = NULL)
{
  check_parms(parms)
  check_string(item_property)
  check_df(domains, c('item_id', item_property))
  check_df(design, 'item_id', nullable=TRUE)
  

  if(length(unique(domains$item_id)) < nrow(domains))
    stop('column domains$item_id must be unique')
  
  if(is.null(design))
  {
    if(is.null(parms$inputs$design))
    {
      # to do: check mst vignettes, if this is not used, better make it an error
      if(inherits(parms,'mst_enorm'))
        message('Computing non-mst profile_tables over an mst design, did you mean to use profile_tables_mst?')
      design = lapply(parms$inputs$bkList, function(bk) tibble(booklet_id=bk$booklet,item_id=bk$items)) %>% 
        bind_rows()
      
      design$item_id = ffactor(design$item_id)
    } else
    {
      design = parms$inputs$design
    }
  }  
  
  domains$item_id = ffactor(as.character(domains$item_id), levels = levels(design$item_id))
  
  if(all(is.na(domains$item_id)))
    stop("none of the item id's in domains occurs in your parameters")
  
  if(anyNA(domains$item_id))
  {
    warning("some item_id's in domains do not occur in your parameters, these have been removed. 
            You can use coef(parms) to see which items are in your parameters")
    domains = filter(domains, !is.na(.data$item_id))
  }

  
  if(!'booklet_id' %in% colnames(design)) 
    design$booklet_id = 'all_items'
  # to do: what if domains not a superset of design? empty property?
  

  profile_tables_(parms, domains, item_property, design) %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
}

profile_tables_ = function(parms, domains, item_property, design)
{
  
  if(!is.factor(domains$item_id))
    domains$item_id = ffactor(as.character(domains$item_id), levels = levels(design$item_id))
  
  design = design %>% 
    select(.data$booklet_id, .data$item_id) %>%
    inner_join(domains[,c('item_id',item_property)], by='item_id') %>%
    mutate(dcat = dense_rank(.data[[item_property]]))
  
  dcat = distinct(design, .data[[item_property]], dcat )
  
  b = if.else(parms$inputs$method == "Bayes", colMeans(parms$est$b), parms$est$b)
  a = parms$inputs$ssIS$item_score

  first = parms$inputs$ssI$first
  last = parms$inputs$ssI$last
  
  parms$inputs$ssI %>% 
    mutate(i = row_number()) %>%
    inner_join(design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    do(
        {
          prof_lmx = E_profile(b,a,first,last, split(.$i, .$dcat))
          
          lapply(prof_lmx, function(mtx)
          {
            tibble(dcat=1:ncol(mtx), expected_domain_score = mtx[1,] )}
          ) %>%
            bind_rows(.id='booklet_score') %>%
            mutate(booklet_score = as.integer(.data$booklet_score)-1L) %>%
            inner_join(dcat, by='dcat')
        }
      ) %>%
      ungroup() %>%
      select(-dcat)
}



# A profile is a table with two rows (earned, not earned) and number of columns
# equal to the number of subsets. As in Norman's Cito brochure.

# Expected profile on subsets A given total score on all items in m
# A is a vector of lists containing mutually exclusive subsets
# subsets defined by the indices of items
# prms is produced by fit_inter
# output is list of profiles for each total score. For score s take E_RM[[s+1]]

# first and last are reduced to only the items occurring in A
# NA categories can be applied by removing them from A as well as first and last
E_profile = function(b, a, first, last, A, cIM=NULL)
{
  nSub = length(A)
  
  items = unlist(A)
  Msc = sum(a[last[items]])
  Msc_sub = integer(nSub)
  
  
  E = matrix(0, nSub, Msc+1)
  for (j in 1:nSub)
  {
    hh = SSTable(b, a, first, last, AB=list(A[[j]], setdiff(items, A[[j]])), cIM=cIM)
      
    Msc_sub[j] = nrow(hh) - 1L
    for (i in 1:nrow(hh))
    {
      for (h in 1:ncol(hh))
      {
        s = i + h - 1L
        sA = i - 1
        E[j,s] = E[j,s] + sA * hh[i,h]
      }
    }   
  }
  Etab = vector("list", Msc+1)
  for (i in 1:length(Etab))
  {
    Etab[[i]] = matrix(0,2,nSub)
    Etab[[i]][1,] = E[,i]
    Etab[[i]][2,] = Msc_sub - E[,i]
  }
  
  Etab
}


# Chi-square disctance between observed and expected
profile_dist <-function(obs, expt)
{
  return(sum((obs-expt)^2/expt))
}
