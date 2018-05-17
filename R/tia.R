

############################
#' Simple test-item analysis
#'
#' Show simple Classical Test Analysis statistics
#' at item and test level
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param type How to present the item level statistics: \code{raw} for each test booklet 
#' separately, \code{averaged} averaged over the test booklet in which the item is included,
#' with the number of persons as weights, or {compared}, in which case the pvalues, 
#' correlations with the sum score (rit), and correlations with the rest score (rit) are 
#' shown in separate tables and compared across booklets
#' @return A list containing:
#' \item{testStats}{a data frame of statistics at test level} 
#' \item{itemStats}{a data frame of statistics at item level}.
#'
tia_tables <- function(dataSrc, predicate = NULL, type=c('raw','averaged','compared')) {
  type = match.arg(type)
  
  qtpredicate = eval(substitute(quote(predicate)))
  x = get_resp_data(dataSrc, qtpredicate, env=caller_env(), summarised=FALSE)$x
  
  # if dataSrc is a db and no predicate then resp_data$design also contains the item_position, might be useful
  
  if(nrow(x) == 0) stop('no data to analyse')
  
  maxScores = x %>%
    group_by(.data$item_id) %>%
    summarise(maxScore=max(.data$item_score)) %>%
    ungroup()

  # suppressed warnings since it will be obvious in the tables from the cor=NA and sd=0
  ti = x %>%
     group_by(.data$booklet_id,.data$item_id) %>%
     summarise(meanScore=mean(.data$item_score),
               sdScore=sd(.data$item_score),
               rit = suppressWarnings(cor(.data$item_score, .data$sumScore)),
               rir = suppressWarnings(cor(.data$item_score, .data$sumScore - .data$item_score)),
               n=n()) %>%
     ungroup() %>%
     inner_join(maxScores, by='item_id') %>%
     mutate(pvalue=.data$meanScore/.data$maxScore)
  

  itemStats = switch(type,
                     raw={
                       ti %>% select(.data$booklet_id, .data$item_id, .data$meanScore, .data$sdScore, .data$maxScore, .data$pvalue, .data$rit, .data$rir, .data$n)
                     },
                     averaged={
                       ti %>% 
                         group_by(.data$item_id) %>%
                         summarise( nBooklets=n(),
                                    meanScore=weighted.mean(.data$meanScore, w=.data$n, na.rm=TRUE),
                                    sdScore=weighted.mean(.data$sdScore, w=.data$n, na.rm=TRUE),
                                    maxScore = max(.data$maxScore),
                                    pvalue=weighted.mean(.data$pvalue, w=.data$n, na.rm=TRUE),
                                    rit=weighted.mean(.data$rit, w=.data$n, na.rm=TRUE),
                                    rir=weighted.mean(.data$rir, w=.data$n, na.rm=TRUE),
                                    n=sum(.data$n, na.rm=TRUE)
                                    ) 
                     },
                     compared={
                       list(
                         pvalue = ti %>% select(.data$booklet_id,.data$item_id,.data$pvalue) %>% 
                                    spread_(key_col='booklet_id',value_col='pvalue'),
                         
                         rit = ti %>% select(.data$booklet_id,.data$item_id,.data$rit) %>% 
                                 spread_(key_col='booklet_id',value_col='rit'),
                         
                         rir = ti %>% select(.data$booklet_id,.data$item_id,.data$rir) %>% 
                                 spread_(key_col='booklet_id',value_col='rir')
                       )
                     }
    )
  testStats = ti %>%
    group_by(.data$booklet_id) %>% 
    filter(complete.cases(.data$sdScore)) %>%
    summarise(nItems=n(),
              alpha=.data$nItems/(.data$nItems-1)*(1-sum(.data$sdScore^2) / sum(.data$rit * .data$sdScore)^2 ),
              meanP=mean(.data$pvalue),
              meanRit=mean(.data$rit),
              meanRir=mean(.data$rir),
              maxTestScore=sum(.data$maxScore),
              N=max(.data$n)) 
  
  list(itemStats=as.data.frame(itemStats), testStats=as.data.frame(testStats))
}

