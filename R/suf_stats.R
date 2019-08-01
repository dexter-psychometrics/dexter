


get_sufStats_nrm = function(respData)
{
  mx = max(respData$x$item_score)
  
  #if this happens to be more than the number of distinct items it does not matter much
  # do not use distinct(item_id), there may be gaps
  nit = length(levels(respData$x$item_id)) 
  
  sufs = suf_stats_nrm(respData$x$booklet_id, respData$x$booklet_score, respData$x$item_id, respData$x$item_score,
                      nit, mx)
  
  # std vectors in c so levels have to be set here
  class(sufs$ssIS$item_id) = 'factor'
  levels(sufs$ssIS$item_id) = levels(respData$x$item_id)
  
  sufs
}


# to do: unit test of this
# # if the score 0 does not occur for an item, it is added with sufI=0 and sufC=0
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
  # to do: mssing 0 score a problem?
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