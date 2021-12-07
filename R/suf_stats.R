


get_sufStats_nrm = function(respData, check_sanity=TRUE)
{
  mx = max(respData$x$item_score)
  if(is.na(mx))
    stop('there is a problem with your response data or database')
  
  design = respData$design
  
  # to~do: I think we need tests with at least three items, fischer criterium
  if(check_sanity)
  {
    if(nrow(design) == 1) 
      stop('There are responses to only one item in your selection, this cannot be calibrated.') 
    
    if(!is_connected(design))
      stop('Your design is not connected')   
  }
  
  #if this happens to be more than the number of distinct items it does not matter much
  # do not use distinct(item_id), there may be gaps
  nit = nlevels(respData$x$item_id)
  
  sufs = suf_stats_nrm(respData$x$booklet_id, respData$x$booklet_score, respData$x$item_id, respData$x$item_score,
                      nit, mx)
  
  # std vectors in c so levels have to be set here
  class(sufs$ssIS$item_id) = 'factor'
  levels(sufs$ssIS$item_id) = levels(respData$x$item_id)
  
  ssIS = sufs$ssIS
  
  # bug in dplyr, min/max of integer in group_by becomes double
  ssI  = ssIS %>% 
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>%
    summarise(first = as.integer(min(.data$rn)),last = as.integer(max(.data$rn))) %>%
    ungroup() %>%
    arrange(.data$item_id)
  
  
  if(check_sanity)
  {
    err = FALSE
    if(any(ssI$first == ssI$last)) 
    {
      message('Items without score variation:')
      print(as.character(ssI$item_id[ssI$first == ssI$last]))
      err = TRUE
    }
    if(any(ssIS$item_score[ssI$first] !=0))
    {

      message('Items without a zero score category')
      print(as.character(ssI$item_id[ssIS$item_score[ssI$first] !=0]))
      err = TRUE
    }
    if(err)
      stop('Some items have only one response category or lack a zero score category',call.=FALSE)
  }
  sufs$ssI = ssI
  
  sufs$design = design %>%
    inner_join(ssI,by='item_id') %>%
    arrange(.data$booklet_id, .data$first)
  
  itm_max = sufs$ssIS %>% 
    group_by(.data$item_id) %>% 
    summarise(maxScore = as.integer(max(.data$item_score))) %>% 
    ungroup()
  
  # max booklet scores
  maxScores = itm_max %>%
    inner_join(sufs$design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    summarise(maxTotScore = sum(.data$maxScore))
  
  # booklets 0:maxscore
  all_scores = maxScores %>% 
    group_by(.data$booklet_id) %>%
    do({tibble(booklet_score=0:.$maxTotScore)}) %>%
    ungroup()
  
  sufs$scoretab = sufs$plt %>%
    distinct(.data$booklet_id, .data$booklet_score,.data$N) %>%
    right_join(all_scores, by=c('booklet_id','booklet_score')) %>%
    mutate(N=coalesce(.data$N, 0L)) %>%
    arrange(.data$booklet_id, .data$booklet_score)
  

  sufs
}



# 
# ssIS = respData$x %>%
#   group_by(.data$item_id, .data$item_score) %>%
#   summarise(sufI=n(), sufC_ = sum(.data$item_score * .data$booklet_score)) %>%
#   ungroup() %>%
#   full_join(tibble(item_id=respData$design$item_id, item_score=0L), by = c("item_id","item_score")) %>%
#   mutate(sufI = coalesce(.data$sufI, 0L), sufC_ = coalesce(.data$sufC_,0L)) %>%
#   arrange(.data$item_id, .data$item_score)
# 
# plt = respData$x %>%
#   group_by(.data$item_id, .data$booklet_score) %>%
#   summarise(meanScore = mean(.data$item_score), N = n()) %>%
#   ungroup()

# change in regard to above implementation is that zero scores are no longer added.
# items without zero score category will now get an error in fit_inter
get_sufStats_im = function(respData)
{
  if(n_distinct(respData$design$booklet_id) != 1)
    stop("invalid resp data for interaction model")
  
  mx = max(respData$x$item_score)
  
  #if this happens to be more than the number of distinct items it does not matter much
  # do not use distinct(item_id), there may be gaps
  nit = length(levels(respData$x$item_id)) 
  
  sufs = suf_stats_im(respData$x$booklet_score, respData$x$item_id, respData$x$item_score,nit, mx)
  
  # std vectors in c so levels have to be set here
  class(sufs$ssIS$item_id) = 'factor'
  levels(sufs$ssIS$item_id) = levels(respData$x$item_id)
  
  class(sufs$plt$item_id) = 'factor'
  levels(sufs$plt$item_id) = levels(respData$x$item_id)

  sufs
}


# ti = respData$x %>%
#   group_by(.data$booklet_id, .data$item_id) %>%
#   summarise(meanScore = mean(.data$item_score),
#             maxScore = max(.data$item_score),
#             sdScore = sd(.data$item_score),
#             rit = suppressWarnings(cor(.data$item_score, .data$booklet_score)),
#             rir = suppressWarnings(cor(.data$item_score, .data$booklet_score - .data$item_score)),
#             n=n()) %>%
#   ungroup() %>%
#   group_by(.data$item_id) %>%
#   mutate(maxScore = max(.data$maxScore)) %>%
#   ungroup() 

get_sufStats_tia = function(respData)
{
  # take length of levels as protection for out of range indexing
  nb = length(levels(respData$design$booklet_id))
  nit = length(levels(respData$design$item_id))
  
  frst_item = respData$design %>%
    distinct(.data$booklet_id, .keep_all=TRUE) %>%
    arrange(.data$booklet_id) %>%
    pull(.data$item_id) %>%
    as.integer()
  
  # indexing in c, make first element empty
  frst_item = c(-10L, frst_item)
  
  tia_C(respData$x$booklet_id, respData$x$booklet_score, respData$x$item_id, respData$x$item_score, 
        nb, nit,
        frst_item, respData$design$booklet_id, respData$design$item_id ) 
}

