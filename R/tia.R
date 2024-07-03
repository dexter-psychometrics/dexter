

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
#' @param distractor add a tia for distractors, only useful for selected response (MC) items
#' @return A list containing:
#' \item{booklets}{a data.frame of statistics at booklet level} 
#' \item{items}{a data.frame (or list if type='compared') of statistics at item level}
#' \item{distractors}{a data.frame of statistics at the response level (if distractor==TRUE), i.e. 
#' rvalue (pvalue for response) and rar (rest-alternative correlation)}
#'
tia_tables = function(dataSrc, predicate = NULL, type=c('raw','averaged','compared'),
                      max_scores = c('observed','theoretical'), distractor=FALSE) 
{
  type = match.arg(type)
  max_scores = match.arg(max_scores)
  check_dataSrc(dataSrc)
  
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  out = list()
  
  respData = get_resp_data(dataSrc, qtpredicate, env=env, summarised=FALSE,
                           extra_columns = if(distractor){'response'}else{NULL})
  
  items = get_sufStats_tia(respData) 
  
  if(max_scores=='theoretical' && is_db(dataSrc))
  {
    items$item_id = as.character(items$item_id)
    items = inner_join(
      select(items, -'max_score'),
      dbGetQuery(dataSrc, 'SELECT item_id, MAX(item_score) AS max_score FROM dxscoring_rules GROUP BY item_id;'),
      by='item_id')
  } 

  items$pvalue = coalesce(items$mean_score / items$max_score, 0)
  
  if(anyNA(items$rit))
    warning("Items without score variation have been removed from the test statistics")
  
  out$booklets = items |>
    filter(complete.cases(.data$rit)) |>
    group_by(.data$booklet_id) |> 
    summarise(n_items=n(),
              alpha=.data$n_items/(.data$n_items-1)*(1-sum(.data$sd_score^2) / sum(.data$rit * .data$sd_score)^2 ),
              mean_pvalue = mean(.data$pvalue),
              mean_rit = mean(.data$rit),
              mean_rir = mean(.data$rir),
              max_booklet_score = sum(.data$max_score),
              n_persons = max(.data$n_persons)) |>
    ungroup() |>
    mutate_if(is.factor, as.character) |>
    df_format()
  
  
  # for presentation purposes, the sd of the item score should be divided by n-1
  # since that is the default in R.
  # Note that this happens AFTER alpha is computed and BEFORE any sd's are grouped over booklets
  
  items$sd_score = sqrt(items$n_persons/(items$n_persons-1)) * items$sd_score
  
  # different views of item statistics
  if(type=='raw')
  {
    out$items = select(items, 'booklet_id', 'item_id', 'mean_score', 'sd_score', 
                           'max_score', 'pvalue', 'rit', 'rir', 'n_persons') |>
      mutate_if(is.factor, as.character) |>
      df_format()
    
  } else if(type=='averaged')
  {
    out$items = items |> 
      group_by(.data$item_id) |>
      summarise( n_booklets = n(),
                 w_mean_score=weighted.mean(.data$mean_score, w = .data$n_persons),
                 sd_score = sqrt(combined_var(.data$mean_score, .data$sd_score^2, .data$n_persons)),
                 max_score = max(.data$max_score),
                 pvalue = weighted.mean(.data$pvalue, w=.data$n_persons),
                 rit = weighted.mean(.data$rit, w=.data$n_persons, na.rm=TRUE),
                 rir = weighted.mean(.data$rir, w=.data$n_persons, na.rm=TRUE),
                 n_persons = sum(.data$n_persons)) |>
      ungroup() |>
      mutate_if(is.factor, as.character) |>
      rename(mean_score = 'w_mean_score') |>
      df_format()
  } else
  {
    items = mutate_if(items, is.factor, as.character)
    
    out$items = list(
      pvalue = items |> 
        select('booklet_id', 'item_id', 'pvalue') |> 
        pivot_wider(names_from='booklet_id', values_from='pvalue', names_sort=TRUE),
      
      rit = items |> 
        select('booklet_id', 'item_id', 'rit') |> 
        pivot_wider(names_from='booklet_id', values_from='rit', names_sort=TRUE),
      
      rir = items |> 
        select('booklet_id', 'item_id', 'rir') |> 
        pivot_wider(names_from='booklet_id', values_from='rir', names_sort=TRUE)
    )
  }
  
  if(distractor)
  {
    d = respData$x |>
      mutate(bs=.data$booklet_score-.data$item_score) |>
      group_by(.data$booklet_id, .data$item_id, .data$response) |>
      summarise(item_score=first(.data$item_score), n=n(), rbsum=sum(.data$bs), rb2sum=sum(.data$bs^2), 
                bsum = sum(.data$booklet_score), b2sum=sum(.data$booklet_score^2), 
                .groups='drop_last') |>
      mutate(N = sum(.data$n), 
             rbmean = sum(.data$rbsum)/.data$N, 
             rb2mean = sum(.data$rb2sum)/.data$N,
             bmean = sum(.data$bsum)/.data$N, 
             b2mean = sum(.data$b2sum)/.data$N,
             rvalue = .data$n/.data$N,
             rar = (.data$rbsum/.data$N - .data$rvalue*.data$rbmean)/
               sqrt(.data$rvalue*(1-.data$rvalue)*(.data$rb2mean - .data$rbmean^2)),
             rat = (.data$bsum/.data$N - .data$rvalue*.data$bmean)/
               sqrt(.data$rvalue*(1-.data$rvalue)*(.data$b2mean - .data$bmean^2))) |>
      ungroup() 
    
    # type==compared makes little sense to me for distractors, so treated same as raw
    if(type=='averaged')
    {
      d = d |>
        group_by(.data$item_id, .data$response, .data$item_score) |>
        summarise(n=sum(.data$n),
                  rvalue = weighted.mean(.data$rvalue,.data$N),
                  rar = weighted.mean(.data$rar,.data$N, na.rm=TRUE),
                  rat = weighted.mean(.data$rat,.data$N, na.rm=TRUE)) |>
        ungroup()
    } else
    {
      d = select(d, 'booklet_id', 'item_id', 'response', 'item_score', 'n', 'rvalue', 'rar','rat')
    }
    
    out$distractors = d |>
      mutate_if(is.factor, as.character) |>
      df_format()
  }
  out
}
