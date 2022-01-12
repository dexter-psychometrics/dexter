

############################
#' Simple test-item analysis
#'
#' Show simple Classical Test Analysis statistics
#' at item and test level
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param type How to present the item level statistics: \code{raw} for each test booklet 
#' separately, \code{averaged} booklets are ignored, with the exception of rit and rir which are averaged over the test booklets,
#' with the number of persons as weights, or \code{compared}, in which case the pvalues, 
#' correlations with the sum score (rit), and correlations with the rest score (rit) are 
#' shown in separate tables and compared across booklets
#' @param max_scores use the observed maximum item score or the theoretical maximum item score 
#' according to the scoring rules in the database to determine pvalues and maximum scores
#' @return A list containing:
#' \item{booklets}{a data.frame of statistics at booklet level} 
#' \item{items}{a data.frame (or list if type='compared') of statistics at item level}
#'
tia_tables = function(dataSrc, predicate = NULL, type=c('raw','averaged','compared'),
                      max_scores = c('observed','theoretical')) 
{
  type = match.arg(type)
  max_scores = match.arg(max_scores)
  check_dataSrc(dataSrc)
  
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, env=env, summarised=FALSE)
  
  items = get_sufStats_tia(respData) 
  
  if(max_scores=='theoretical' && is_db(dataSrc))
  {
    items$item_id = as.character(items$item_id)
    items = inner_join(
      select(items, -.data$max_score),
      dbGetQuery(dataSrc, 'SELECT item_id, MAX(item_score) AS max_score 
                          FROM dxscoring_rules GROUP BY item_id;'),
      by='item_id')
  } 

  items$pvalue = coalesce(items$mean_score / items$max_score, 0)
  
  if(anyNA(items$rit))
    warning("Items without score variation have been removed from the test statistics")
  
  booklets = items %>%
    filter(complete.cases(.data$rit)) %>%
    group_by(.data$booklet_id) %>% 
    summarise(n_items=n(),
              alpha=.data$n_items/(.data$n_items-1)*(1-sum(.data$sd_score^2) / sum(.data$rit * .data$sd_score)^2 ),
              mean_pvalue = mean(.data$pvalue),
              mean_rit = mean(.data$rit),
              mean_rir = mean(.data$rir),
              max_booklet_score = sum(.data$max_score),
              n_persons = max(.data$n_persons)) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
  
  # different views of item statistics
  if(type=='raw')
  {
    items = select(items, .data$booklet_id, .data$item_id, .data$mean_score, .data$sd_score, 
                           .data$max_score, .data$pvalue, .data$rit, .data$rir, .data$n_persons) %>%
      mutate_if(is.factor, as.character) %>%
      df_format()
    
  } else if(type=='averaged')
  {
    items = items %>% 
      group_by(.data$item_id) %>%
      summarise( n_booklets = n(),
                 w_mean_score=weighted.mean(.data$mean_score, w = .data$n_persons),
                 sd_score = sqrt(combined_var(.data$mean_score, .data$sd_score^2, .data$n_persons)),
                 max_score = max(.data$max_score),
                 pvalue = weighted.mean(.data$pvalue, w=.data$n_persons),
                 rit = weighted.mean(.data$rit, w=.data$n_persons, na.rm=TRUE),
                 rir = weighted.mean(.data$rir, w=.data$n_persons, na.rm=TRUE),
                 n_persons = sum(.data$n_persons)) %>%
      ungroup() %>%
      mutate_if(is.factor, as.character) %>%
      rename(mean_score=.data$w_mean_score) %>%
      df_format()
  } else
  {
    items = mutate_if(items, is.factor, as.character)
    
    items = list(
      pvalue = items %>% 
        select(.data$booklet_id,.data$item_id,.data$pvalue) %>% 
        spread(key=.data$booklet_id, value=.data$pvalue),
      
      rit = items %>% 
        select(.data$booklet_id,.data$item_id,.data$rit) %>% 
        spread(key=.data$booklet_id,value=.data$rit),
      
      rir = items %>% 
        select(.data$booklet_id,.data$item_id,.data$rir) %>% 
        spread(key=.data$booklet_id,value=.data$rir)
    )
  }
  
  
  list(booklets=booklets, items=items)
}

