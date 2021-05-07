

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
#' separately, \code{averaged} averaged over the test booklet in which the item is included,
#' with the number of persons as weights, or \code{compared}, in which case the pvalues, 
#' correlations with the sum score (rit), and correlations with the rest score (rit) are 
#' shown in separate tables and compared across booklets
#' @param max_scores use the observed maximum item score or the theoretical maximum items score 
#' according to the scoring rules in the database to compute pvalues and maximum scores
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
  env=caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, env=env, summarised=FALSE)
  
  ti = get_sufStats_tia(respData) 
  
  if(max_scores=='theoretical' && is_db(dataSrc))
  {
    ti$item_id = as.character(ti$item_id)
    ti = inner_join(
      select(ti,-.data$maxScore),
      dbGetQuery(dataSrc, 'SELECT item_id, MAX(item_score) AS max_score 
                          FROM dxscoring_rules GROUP BY item_id;'),
      by='item_id') %>%
      rename(maxScore = .data$max_score)
  } 
  
  ti = mutate(ti, pvalue=coalesce(.data$meanScore/.data$maxScore, 0))
  
  
  if(type=='raw')
  {
    itemStats = select(ti, .data$booklet_id, .data$item_id, .data$meanScore, .data$sdScore, 
                           .data$maxScore, .data$pvalue, .data$rit, .data$rir, .data$n) %>%
      mutate_if(is.factor, as.character) %>%
      df_format()
    
  } else if(type=='averaged')
  {
    itemStats = ti %>% 
      group_by(.data$item_id) %>%
      summarise( nBooklets=n(),
                 meanScore=weighted.mean(.data$meanScore, w=.data$n),
                 sdScore=weighted.mean(.data$sdScore, w=.data$n),
                 maxScore = max(.data$maxScore),
                 pvalue=weighted.mean(.data$pvalue, w=.data$n),
                 rit=weighted.mean(.data$rit, w=.data$n, na.rm=TRUE),
                 rir=weighted.mean(.data$rir, w=.data$n, na.rm=TRUE),
                 n=sum(.data$n)) %>%
      ungroup() %>%
      mutate_if(is.factor, as.character) %>%
      df_format()
  } else
  {
    ti = mutate_if(ti, is.factor, as.character)
    itemStats = list(
      pvalue = ti %>% 
        select(.data$booklet_id,.data$item_id,.data$pvalue) %>% 
        spread(key=.data$booklet_id, value=.data$pvalue),
      
      rit = ti %>% 
        select(.data$booklet_id,.data$item_id,.data$rit) %>% 
        spread(key=.data$booklet_id,value=.data$rit),
      
      rir = ti %>% 
        select(.data$booklet_id,.data$item_id,.data$rir) %>% 
        spread(key=.data$booklet_id,value=.data$rir)
    )
  }
  if(anyNA(ti$rit))
    warning("Items without score variation have been removed from the test statistics")

  testStats = ti %>%
    filter(complete.cases(.data$rit)) %>%
    group_by(.data$booklet_id) %>% 
    summarise(nItems=n(),
              alpha=.data$nItems/(.data$nItems-1)*(1-sum(.data$sdScore^2) / sum(.data$rit * .data$sdScore)^2 ),
              meanP=mean(.data$pvalue),
              meanRit=mean(.data$rit),
              meanRir=mean(.data$rir),
              maxTestScore=sum(.data$maxScore),
              N=max(.data$n)) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
  
  # new names since version 1.0.3 may 2021, old are kept for backward compatibility
  booklets = rename(testStats, n_items=.data$nItems, mean_pvalue=.data$meanP,
                    mean_rit=.data$meanRit, mean_rir=.data$meanRir, 
                    max_booklet_score=.data$maxTestScore, n_persons=.data$N)
  
  items = rename(itemStats, mean_score = .data$meanScore,
                 sd_score=.data$sdScore, max_score=.data$maxScore, n_persons=.data$n)
  
  list(booklets=booklets,items=items,itemStats=itemStats, testStats=testStats)
}

