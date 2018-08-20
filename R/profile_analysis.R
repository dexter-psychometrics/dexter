


##########################################
#' Profile analysis
#'
#' Expected and observed domain scores, conditional on the test score, per person or test score. 
#' Domains are specified as categories of items using item_properties.
#'
#' @param dataSrc a dexter project db handle or a data.frame with columns: person_id, item_id, item_score, 
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
#' \item{profiles}{a data.frame with columns person_id, booklet_id, sumScore, 
#' -item_property-, domain_score, expected_domain_score}
#' \item{profile_tables}{a data.frame with columns booklet_id, sumScore, 
#' -item_property-, expected_domain_score }
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
profiles = function(dataSrc, parms, item_property, predicate=NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  
  
  if(!inherits(parms,'prms'))
    stop('argument `parms` should be a parameters object resulting from fit_enorm')
  
  if(inherits(dataSrc, 'data.frame'))
  {
    if(!item_property %in% colnames(dataSrc))
      stop(paste0('dataSrc should include a column `',item_property,'`'))
    
    idom = distinct(dataSrc, .data$item_id, .data[[!item_property]])
    
    if(nrow(idom) > nrow(distinct(idom, .data$item_id)))
      stop('Each item may only have a single item property')
  } else
  {
    item_property = tolower(item_property)
    
    if(!item_property %in% dbListFields(dataSrc,'dxitems'))
      stop(paste0('There is no item_property with name `',item_property,'` in your project'))
  }
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=caller_env(), 
                           extra_columns = item_property, extra_design_columns = item_property)
  
  respData$x %>%
    group_by(.data$person_id, .data$booklet_id, .data$sumScore, .data[[!!item_property]]) %>%
    summarise(domain_score = sum(.data$item_score)) %>%
    ungroup() %>%
    inner_join(
        profile_tables(parms = parms, design = respData$design, 
                       domains=distinct(respData$design, .data$item_id, .data[[!!item_property]]),
                       item_property = item_property),
        by = c('booklet_id','sumScore',item_property)) %>%
    as.data.frame()

}
  
#' @rdname profiles
profile_tables = function(parms, domains, item_property, design = NULL)
{
  if(!inherits(parms,'prms'))
    stop('argument `parms` should be a parameters object resulting from fit_enorm')
  if(!inherits(domains,'data.frame'))
    stop('argument domains should be a data.frame')
  if(!'item_id' %in% colnames(domains))
    stop('domains should include a column `item_id`')
  if(!item_property %in% colnames(domains))
    stop(paste0('domains should include a column `',item_property,'`'))
  if(length(unique(domains$item_id)) < nrow(domains))
    stop('column domains$item_id must be unique')
  
  
  if(is.null(design))
  {
    if(is.null(parms$inputs$design))
    {
      if(inherits(parms,'mst_enorm'))
        message('Computing non-mst profile_tables over an mst design, did you mean to use profile_tables_mst?')
      design = lapply(parms$inputs$bkList, function(bk) tibble(booklet_id=bk$booklet,item_id=bk$items)) %>% bind_rows()
    } else
    {
      design = parms$inputs$design
    }
  }  

  
  if(!'booklet_id' %in% colnames(design)) design$booklet_id = 'all_items'
  # what if domains not a superset of design?
  
  
  design = design %>% 
    select(.data$booklet_id, .data$item_id) %>%
    inner_join(domains, by='item_id') %>%
    mutate(dcat = dense_rank(.data[[!!item_property]]))
  
  dcat = distinct(design, .data[[!!item_property]], dcat )
  
    parms$inputs$ssI %>% 
      mutate(i = row_number()) %>%
      inner_join(design, by='item_id') %>%
      group_by(.data$booklet_id) %>%
      do(
        {
          prof_lmx = E_profile_enorm(parms, split(.$i, .$dcat))
          
          lapply(prof_lmx$E_RM, function(mtx)
          {
            tibble(dcat=1:ncol(mtx), expected_domain_score = mtx[1,] )}
          ) %>%
            bind_rows(.id='sumScore') %>%
            mutate(sumScore = as.integer(.data$sumScore)-1L) %>%
            inner_join(dcat, by='dcat')
        }
      ) %>%
      ungroup() %>%
      select(-dcat) %>%
      as.data.frame()

}

# prms is enorm params
# A is a vector of lists containing mutually exclusive subsets (indexes in ssI)
# TO DO; consider possibilities of Bayesian parms
E_profile_enorm <- function(prms, A)
{
  nSub=length(A)
  if (prms$inputs$method == "Bayes") prms$est$b = colMeans(prms$est$b)
  Msc = sum(prms$inputs$ssIS$item_score[prms$inputs$ssI[sort(unlist(A)),]$last])
  Msc_sub = rep(0,nSub)
  
  AA=vector("list",2)
  E_RM = matrix(0, nSub, Msc+1)
  for (j in 1:nSub)
  {
    AA[[1]] = A[[j]]
    AA[[2]] = unlist(A[setdiff(1:nSub,j)])
    hh_RM = SSTable_ENORM(prms$est$b, prms$inputs$ssIS$item_score, prms$inputs$ssI$first, prms$inputs$ssI$last, AA)$tbl

    Msc_sub[j] = nrow(hh_RM)-1
    for (i in 1:nrow(hh_RM))
    {
      for (h in 1:ncol(hh_RM))
      {
        s = i+h-1
        sA = i-1
        E_RM[j,s] = E_RM[j,s] + sA*hh_RM[i,h]
      }
    }
  }
  Etab_RM = vector("list", Msc+1)
  for (i in 1:length(Etab_RM))
  {
    Etab_RM[[i]] = matrix(0,2,nSub)
    Etab_RM[[i]][1,] = E_RM[,i]
    Etab_RM[[i]][2,] = Msc_sub - E_RM[,i]
  }
  return(list(E_RM=Etab_RM))
}


  
# A profile is a table with two rows (earned, not earned) and number of columns
# equal to the number of subsets. As in Norman's Cito brochure.

# Expected profile on subsets A given total score on all items in m
# A is a vector of lists containing mutually exclusive subsets
# subsets defined by the indices of items
# prms is produced by fit_inter
# output is list of profiles for each total score. For score s take E_RM[[s+1]]
E_profile <- function(prms, A)
{
  nSub=length(A)
  Msc = sum(prms$ss$sl$item_score[prms$ss$il$last])
  Msc_sub = rep(0,nSub)
  
  AA=vector("list",2)
  E_RM = matrix(0, nSub, Msc+1)
  E_IM = matrix(0, nSub, Msc+1)
  for (j in 1:nSub)
  {
    AA[[1]] = A[[j]]
    AA[[2]] = unlist(A[setdiff(1:nSub,j)])
    hh_RM = SSTable(prms,AA,"RM")$tbl
    hh_IM = SSTable(prms,AA,"IM")$tbl
    Msc_sub[j] = nrow(hh_RM)-1
    for (i in 1:nrow(hh_RM))
    {
      for (h in 1:ncol(hh_RM))
      {
        s=i+h-1
        sA=i-1
        E_RM[j,s] = E_RM[j,s] + sA*hh_RM[i,h]
        E_IM[j,s] = E_IM[j,s] + sA*hh_IM[i,h] 
      }
    }
  }
  Etab_RM = vector("list", Msc+1)
  Etab_IM = vector("list", Msc+1)
  for (i in 1:length(Etab_RM))
  {
    Etab_RM[[i]] = matrix(0,2,nSub)
    Etab_RM[[i]][1,] = E_RM[,i]
    Etab_RM[[i]][2,] = Msc_sub - E_RM[,i]
    Etab_IM[[i]] = matrix(0,2,nSub)
    Etab_IM[[i]][1,] = E_IM[,i]
    Etab_IM[[i]][2,] = Msc_sub - E_IM[,i]
  }
  return(list(E_RM = Etab_RM, E_IM = Etab_IM))
}

# Chi-square disctance between observed and expected
profile_dist <-function(obs, expt)
{
  return(sum((obs-expt)^2/expt))
}

r_profile <- function(b, a, c, A, first, last, score, model)
{
  x_scores = rep(0,length(A))
  Ei=1:length(first)
  s=score
  it=1
  if (model=="IM")
  {
    logb=log(b)
    logc=log(c)
    while ((s>0)&(it<=length(first)))
    {
      eta = exp(logb + (a * s) * logc)
      g = elsym(eta, a, first[Ei], last[Ei])
      gi = elsym(eta, a, first[Ei[-1]], last[Ei[-1]])
      pi = rep(0,last[it]-first[it]+1)
      k=1
      for (j in first[it]:last[it]) 
      {
        idx = s + 1 - a[j]
        if ((idx > 0) & (idx <= length(gi))) 
        {
          pi[k] = exp(log(eta[j]) + log(gi[idx]) - log(g[s + 1]))
        }
        k=k+1
      }
      resp = sample(a[first[it]:last[it]], 1 , prob=pi)
      which_set=which(unlist(lapply(A,function(x) it%in%x)))
      x_scores[which_set]=x_scores[which_set]+resp
      s = s - a[first[it]:last[it]][resp+1]
      it = it + 1
      Ei = Ei[-1]
    }
  }else
  {
    while ((s>0)&(it<=length(first)))
    {
      g = elsym(b, a, first[Ei], last[Ei])
      gi = elsym(b, a, first[Ei[-1]], last[Ei[-1]])
      pi = rep(0,last[it]-first[it]+1)
      k=1
      for (j in first[it]:last[it]) 
      {
        idx = s + 1 - a[j]
        if ((idx > 0) & (idx <= length(gi))) 
        {
          pi[k] = exp(log(b[j]) + log(gi[idx]) - log(g[s + 1]))
        }
        k=k+1
      }
      resp = sample(a[first[it]:last[it]], 1 , prob=pi)
      which_set=which(unlist(lapply(A,function(x) it%in%x)))
      x_scores[which_set]=x_scores[which_set]+resp
      s = s - a[first[it]:last[it]][resp+1]
      it = it + 1
      Ei = Ei[-1]
    }
  }
  
  ## make profile
  Msc_sub = unlist(lapply(A,function(x) sum(a[last[x]])))
  sim_tab = matrix(0,2,length(A))
  sim_tab[1,] = x_scores
  sim_tab[2,] = Msc_sub - x_scores
  
  return(sim_tab)
}

# For all profiles defined by A, the probabilities given total-score
# under ENORM. Currently only CML
enum_prof <-function(prms, A, sums=FALSE)
{
  nA = length(A)
  Msc = rep(0,nA)
  a = prms$inputs$ssIS$item_score
  b = prms$est$b
  first = prms$inputs$ssI$first
  last = prms$inputs$ssI$last
  for (i in 1:length(A)) Msc[i] = sum(a[last[A[[i]]]])
  indx = unlist(A, use.names = FALSE)
  lg = log(elsym(b, a, first[indx], last[indx]))
  P = array(0,Msc+1)
  SM = P
  for (i in 1:length(A))
  {
    prefix <- paste(rep(",", i - 1), collapse = "")
    suffix <- paste(rep(",", nA - i), collapse = "")
    lgA = log(elsym(b, a, first[A[[i]]], last[A[[i]]]))
    for (j in 0:Msc[i])
    {
      expr <- paste("P[", paste(prefix, j+1, suffix, sep = ""), "] <-",
                    "P[", paste(prefix, j+1, suffix, sep = ""), "]+lgA[j+1]", sep = "")
      eval(parse(text = expr))
      expr <- paste("SM[", paste(prefix, j+1, suffix, sep = ""), "] <-",
                    "SM[", paste(prefix, j+1, suffix, sep = ""), "]+j", sep = "")
      eval(parse(text = expr))
    }
  }
  for (i in 0:sum(Msc))
  {
    indx = which(SM==i, arr.ind = T, useNames = F)
    P[indx] = P[indx] - lg[i+1] 
  }
  if (!sums){
    return(exp(P))
  }else
  {
    return(list(P=exp(P), SM=SM))
  }
}


