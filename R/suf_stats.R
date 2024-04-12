
possible_scores = function(a, first, last) drop(possible_scores_C(as.integer(a), as.integer(first-1L), as.integer(last-1L)))

get_sufStats_nrm = function(respData, check_sanity=TRUE)
{
    mx = max(respData$x$item_score)
    if(is.na(mx))
      stop_('there is a problem with your response data or database')
    
    design = respData$design
    
    # to~do: I think we need tests with at least three items, fischer criterium
    if(check_sanity)
    {
      if(nrow(design) == 1) 
        stop_('There are responses to only one item in your selection, this cannot be calibrated.') 
      
      if(!is_connected(design))
        stop_('Your design is not connected')   
    }
    
    #if this happens to be more than the number of distinct items it does not matter much
    # do not use distinct(item_id), there may be gaps
    nit = nlevels(respData$x$item_id)
    
    sufs = suf_stats_nrm(respData$x$booklet_id, respData$x$booklet_score, respData$x$item_id, respData$x$item_score,
                         nit, mx)
    
    # std vectors in c so levels have to be set here
    class(sufs$ssIS$item_id) = 'factor'
    levels(sufs$ssIS$item_id) = levels(respData$x$item_id)
    
    #sum is necessary oin case 0 cat is missing
    n_all = sufs$ssIS |>
      group_by(.data$item_id) |>
      summarise(n_rsp = sum(.data$sufI), n0=sum(.data$sufI[.data$item_score==0])) |>
      ungroup()
    
    if(check_sanity)
    {
      err = FALSE
      if(any(n_all$n0==0))
      {
        message('Items without a zero score category')
        print(as.character(n_all$item_id[n_all$n0==0]))
        err = TRUE
      }
      if(any(n_all$n_rsp==n_all$n0))
      {
        message('Items without score variation')
        print(as.character(n_all$item_id[n_all$n_rsp==n_all$n0]))
        err = TRUE
      }
      if(err)
        stop_('Some items have only one response category or lack a zero score category')
    }
    
    
    sufs$ssIS = filter(sufs$ssIS,.data$item_score>0)
    
    # bug in dplyr, min/max of integer in group_by becomes double
    sufs$ssI  = sufs$ssIS |> 
      mutate(rn = row_number()) |>
      group_by(.data$item_id) |>
      summarise(first = as.integer(min(.data$rn)),last = as.integer(max(.data$rn))) |>
      ungroup() |>
      inner_join(n_all,by='item_id') |>
      arrange(.data$item_id)
    
    
    sufs$design = design |>
      inner_join(sufs$ssI,by='item_id') |>
      arrange(.data$booklet_id, .data$first)
    
    itm_max = sufs$ssIS |> 
      group_by(.data$item_id) |> 
      summarise(maxScore = as.integer(max(.data$item_score))) |> 
      ungroup()
    
    # max booklet scores
    maxScores = itm_max |>
      inner_join(sufs$design, by='item_id') |>
      group_by(.data$booklet_id) |>
      summarise(maxTotScore = sum(.data$maxScore))
    
    # booklets 0:maxscore
    all_scores = maxScores |> 
      group_by(.data$booklet_id) |>
      do({tibble(booklet_score=0:.$maxTotScore)}) |>
      ungroup()
    
    sufs$scoretab = sufs$plt |>
      distinct(.data$booklet_id, .data$booklet_score,.data$N) |>
      right_join(all_scores, by=c('booklet_id','booklet_score')) |>
      mutate(N=coalesce(.data$N, 0L)) |>
      arrange(.data$booklet_id, .data$booklet_score)
    
    
    m = sufs$scoretab |>
      group_by(.data$booklet_id) |>
      summarise(M=sum(.data$N))
    
    
    sufs$booklet = sufs$design |>
      group_by(.data$booklet_id) |>
      summarise(max_score = sum(sufs$ssIS$item_score[.data$last]),
                nit=n()) |>
      ungroup() |>
      inner_join(m,by='booklet_id') |>
      mutate(n_scores=.data$max_score+1L) |>
      arrange(.data$booklet_id)
    
    sufs
}




# 
# ssIS = respData$x |>
#   group_by(.data$item_id, .data$item_score) |>
#   summarise(sufI=n(), sufC_ = sum(.data$item_score * .data$booklet_score)) |>
#   ungroup() |>
#   full_join(tibble(item_id=respData$design$item_id, item_score=0L), by = c("item_id","item_score")) |>
#   mutate(sufI = coalesce(.data$sufI, 0L), sufC_ = coalesce(.data$sufC_,0L)) |>
#   arrange(.data$item_id, .data$item_score)
# 
# plt = respData$x |>
#   group_by(.data$item_id, .data$booklet_score) |>
#   summarise(meanScore = mean(.data$item_score), N = n()) |>
#   ungroup()

# change in regard to above implementation is that zero scores are no longer added.
# items without zero score category will now get an error in fit_inter
get_sufStats_im = function(respData, check_sanity=TRUE)
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

  sufs$ssI = sufs$ssIS |>
    group_by(.data$item_id) |>
    summarise(nCat = n()-1L, sufC = sum(.data$sufC_), item_maxscore = max(.data$item_score), 
              item_minscore = min(.data$item_score), n_rsp=sum(.data$sufI), n0=sum(.data$sufI[.data$item_score==0])) |>
    ungroup() |>
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1L, last = cumsum(.data$nCat))  |>
    arrange(.data$item_id)
  
 
  
  # theoretical max score on the test
  maxTestScore = sum(sufs$ssI$item_maxscore)
  
  # scoretab, include unachieved and impossible scores
  sufs$scoretab = sufs$plt |>
    select('booklet_score', 'N') |>
    distinct(.data$booklet_score, .keep_all=TRUE) |>
    right_join(tibble(booklet_score=0L:maxTestScore), by="booklet_score") |>
    mutate(N=coalesce(.data$N, 0L)) |>
    arrange(.data$booklet_score)
  
  if(check_sanity)
  {
    if(any(sufs$ssI$nCat<1))
    {
      message('The following items have no score variation:')
      sufs$ssI |>
        filter(.data$nCat<2) |>
        pull(.data$item_id) |>
        as.character() |>
        print()
      stop("data contains items without score variation")
    }
    
    if(any(sufs$ssI$item_minscore>0))
    {
      message('The following items have no zero score category:')
      sufs$ssI |>
        filter(.data$item_minscore>0) |>
        pull(.data$item_id) |>
        as.character() |>
        print()
      stop("data contains items without zero score category")
    }
    if(all_trivial_scores(sufs$ssIS))
      warning("every score can be reached in only one way, no data reduction possible")
  }
  
  sufs$ssIS = filter(sufs$ssIS,.data$item_score>0)
  
  sufs
}









# ti = respData$x |>
#   group_by(.data$booklet_id, .data$item_id) |>
#   summarise(meanScore = mean(.data$item_score),
#             maxScore = max(.data$item_score),
#             sdScore = sd(.data$item_score),
#             rit = suppressWarnings(cor(.data$item_score, .data$booklet_score)),
#             rir = suppressWarnings(cor(.data$item_score, .data$booklet_score - .data$item_score)),
#             n=n()) |>
#   ungroup() |>
#   group_by(.data$item_id) |>
#   mutate(maxScore = max(.data$maxScore)) |>
#   ungroup() 


get_sufStats_tia = function(respData)
{
  # take length of levels as protection for out of range indexing
  nb = length(levels(respData$design$booklet_id))
  nit = length(levels(respData$design$item_id))
  
  frst_item = respData$design |>
    distinct(.data$booklet_id, .keep_all=TRUE) |>
    arrange(.data$booklet_id) |>
    pull(.data$item_id) |>
    as.integer()
  
  # indexing in c, make first element empty
  frst_item = c(-10L, frst_item)
  
  tia_C(respData$x$booklet_id, respData$x$booklet_score, respData$x$item_id, respData$x$item_score, 
        nb, nit,
        frst_item, respData$design$booklet_id, respData$design$item_id ) 

}
