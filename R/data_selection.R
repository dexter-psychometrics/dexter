# TO DO: rlang has an is_reference method, don't have to use home baked one in cpp
# test and use


#' Selecting data
#' 
#' Extract data from a dexter database 
#' 
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate an expression to select data on
#' @param columns the columns you wish to select, can include any column in the project, see: \code{\link{get_variables}}
#' @return a data.frame of responses
#' @details 
#' Many functions in Dexter accept a data source and a predicate. Predicates are extremely flexible 
#' but they have a few limitations because they work on the individual response level. It is therefore not possible
#' for example, to remove complete person cases from an analysis based on responses to a single item 
#' by using just a predicate expression.
#' 
#' For such cases, Dexter supports selecting the data and manipulating it before passing it back to a Dexter function 
#' or possibly doing something else with it. The following example will hopefully clarify this.
#' @examples 
#' \donttest{
#' \dontrun{
#' # goal: fit the extended nominal response model using only persons 
#' # without any missing responses
#' library(dplyr)
#' 
#' # the following would not work since it will omit only the missing 
#' # responses, not the persons; which is not what we want in this case
#' wrong = fit_enorm(db, response != 'NA')
#' 
#' # to select on an aggregate level, we need to gather the data and 
#' # manipulate it ourselves
#' data = get_responses(db, 
#'    columns=c('person_id','item_id','item_score','response')) |>
#'    group_by(person_id) |>
#'    mutate(any_missing = any(response=='NA')) |>
#'    filter(!any_missing)
#'
#' correct = fit_enorm(data)
#' 
#' }}
get_responses = function(dataSrc, predicate=NULL, columns=c('person_id','item_id','item_score'))
{
  env = caller_env() 
  qtpredicate = eval(substitute(quote(predicate)))
  check_dataSrc(dataSrc)
  df_format(get_responses_(dataSrc, qtpredicate=qtpredicate, env=env, columns=columns))
}

get_responses_ = function(dataSrc, qtpredicate=NULL, env=NULL, columns=c('person_id','item_id','item_score'))
{
  if(is.null(env))
    env = caller_env() 

  if(inherits(dataSrc,'data.frame'))
  {
    if(is.null(qtpredicate))
      return(dataSrc[,columns])    
    
    dataSrc[eval_tidy(qtpredicate, data=dataSrc, env=env), columns]
  
  } else if(inherits(dataSrc,'matrix'))
  {
    items = if.else(is.null(colnames(dataSrc)), sprintf('item%04i',1:ncol(dataSrc)), colnames(dataSrc))
    persons = if.else(is.null(rownames(dataSrc)), 1:nrow(dataSrc), rownames(dataSrc))
    
    if(!is.null(qtpredicate))
      stop('predicates not supported for matrix datasource')
    
    tibble(person_id=rep(persons, ncol(dataSrc)), 
                  item_id=rep(items, each=nrow(dataSrc)),
                  item_score=as.integer(dataSrc)) |>
      filter(!is.na(.data$item_score))
    
  } else
  {
    db_get_responses(dataSrc, qtpredicate=qtpredicate, env=env, columns=columns)
  }
  
}


get_cte = function(db, columns)
{
  columns = setdiff(columns, 
    c(dbListFields(db,'dxresponses'), dbListFields(db,'dxscoring_rules')))
  
  cte = " dxresponses INNER JOIN dxscoring_rules USING(item_id, response)"
  
  if(length(intersect(dbListFields(db,'dxpersons'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxpersons USING(person_id)')

  if(length(intersect(dbListFields(db,'dxadministrations'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxadministrations USING(person_id, booklet_id)')
    
  if(length(intersect(dbListFields(db,'dxitems'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxitems USING(item_id)')
  
  if(length(intersect(dbListFields(db,'dxbooklets'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxbooklets USING(booklet_id)')
  
  if(length(intersect(dbListFields(db,'dxbooklet_design'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxbooklet_design USING(booklet_id, item_id)')
  
  paste(cte, collapse=' ')
}

db_get_responses = function(db, qtpredicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) 
    env=caller_env() 
  
  columns = dbValid_colnames(columns)
	pred_sql = qtpredicate_to_sql(qtpredicate, db, env)
  used_columns = union(pred_sql$db_vars, columns)
  
  # decide which tables we need to join since joining costs time
  cte = get_cte(db, used_columns)
  
  if(pred_sql$success)
  {   
    respData = try(dbGetQuery(db, 
                     paste("SELECT", paste0(columns, collapse=','),
                           "FROM", cte, pred_sql$where)),
                   silent=TRUE)
    if(!inherits(respData,'try-error'))
      return(respData)
    
    if(grepl('malformed', respData))
      stop('your database file appears to be broken, this can happen due to a copy error. For more information,', 
           'see: https://www.sqlite.org/faq.html#q21')
  }
  #message('sql translation of the predicate failed')
  #print(pred_sql)
  #print(qtpredicate)
  
  # translation to sql did not work  

  respData = dbGetQuery(db, 
               paste("SELECT", paste(used_columns, collapse=','),
                      "FROM", cte))
  
  qtpredicate = correct_symbol_case(qtpredicate, used_columns,env=env)
  
  return(respData[eval_tidy(qtpredicate, data=respData, env=env), columns])
}

# be careful not to call with unsafe queries
# be careful not to call including columns that are not person or booklet properties
#   as these will be used in the grouping
# if columns includes booklet_id, will not merge over persons
# if columns does not include booklet_id, will merge over persons
# an extra column booklet_score will be returned always
db_get_testscores = function(db, qtpredicate=NULL, columns=c('booklet_id','person_id'), env=NULL)
{
  if(is.null(env)) 
    env=caller_env() 
  
  columns = dbValid_colnames(columns)
  pred_sql = qtpredicate_to_sql(qtpredicate, db, env)
  used_columns = union(pred_sql$db_vars, columns)
  
  # decide which tables we need to join since joining costs time
  cte = get_cte(db, used_columns)
  
  if(pred_sql$success)
  { 
    respData = try(dbGetQuery(db, 
                     paste("SELECT", paste0(columns, collapse=','), ", SUM(item_score) AS booklet_score",
                             "FROM", cte, pred_sql$where,
                             "GROUP BY", paste0(columns, collapse=','),";")),        
                   silent=TRUE) 
    if(!inherits(respData,'try-error'))
      return(respData)
  }
  
  respData = dbGetQuery(db, 
               paste("SELECT", paste(union(used_columns,'item_score'), collapse=','),
                        "FROM", cte))
  
  qtpredicate = correct_symbol_case(qtpredicate, used_columns,env=env)
  
  # worst case scenario
  respData[eval_tidy(qtpredicate, data=respData, env=env), union(columns,'item_score')] |>
      group_by_at(columns) |>
      summarise(booklet_score = sum(.data$item_score)) |>
      ungroup()
}





# returns design with column n_persons added
# optimized for different sorts of predicates
db_get_design = function(db, qtpredicate=NULL, env=NULL)
{
  if(is.null(env))
    env = caller_env()
  
  if(is.null(qtpredicate))
    return(
      dbGetQuery(db,
        "WITH bkl_count AS(
          SELECT booklet_id, COUNT(*) AS n_persons 
            FROM dxadministrations
              GROUP BY booklet_id)
            
         SELECT booklet_id, item_id, item_position, n_persons 
           FROM dxbooklet_design
             INNER JOIN bkl_count USING(booklet_id);"))
  
  if(is_bkl_safe(db, qtpredicate, env) )
  {
    uns = setdiff(c(dbListFields(db,'dxadministrations'), dbListFields(db,'dxpersons')),
                  'booklet_id')
    pred_sql = qtpredicate_to_sql(qtpredicate, db, env)

    if(pred_sql$success)
    {
      if(length(intersect(pred_sql$db_vars, uns)) == 0) 
      {
        res = try(
          dbGetQuery(db, paste(
            "WITH bkl_count AS(
              SELECT booklet_id, COUNT(*) AS n_persons 
                FROM dxadministrations
                  GROUP BY booklet_id)
                
             SELECT booklet_id, item_id, item_position, n_persons 
               FROM dxbooklets
                 INNER JOIN bkl_count USING(booklet_id)
                   INNER JOIN dxbooklet_design USING(booklet_id)
                     INNER JOIN dxitems USING(item_id)",
              pred_sql$where,";")))
      } else
      {
        res = try(
          dbGetQuery(db, paste(
            "SELECT booklet_id, item_id, item_position, COUNT(*) AS n_persons 
               FROM dxbooklets
                 INNER JOIN dxbooklet_design USING(booklet_id)
                   INNER JOIN dxitems USING(item_id)
                    INNER JOIN dxadministrations USING(booklet_id)
                      INNER JOIN dxpersons USING(person_id)",
                pred_sql$where,
            "GROUP BY booklet_id, item_id, item_position;")))
      }
      if(!inherits(res,'try-error'))
        return(res)
    }
  }
  
  # fail safe + case when !booklet_safe, will be slower  
  get_resp_data(db, qtpredicate, env=env, extra_columns = 'item_position')$x |>
    count(.data$booklet_id, .data$item_id, .data$item_position, name='n_persons')

}
