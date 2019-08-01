# to do: postgres apparently has lower case table names  or case sensitive table names

#' Selecting data
#' 
#' Extract data from a dexter database 
#' 
#' @param dataSrc a dexter project database or data.frame
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
#'    columns=c('person_id','item_id','item_score','response')) %>%
#'    group_by(person_id) %>%
#'    mutate(any_missing = any(response=='NA')) %>%
#'    filter(!any_missing)
#'
#' correct = fit_enorm(data)
#' 
#' }}
get_responses = function(dataSrc, predicate=NULL, columns=c('person_id','item_id','item_score'))
{
  env = caller_env() 
  qtpredicate = eval(substitute(quote(predicate)))
  get_responses_(dataSrc, qtpredicate=qtpredicate, env=env, columns=columns)
}

get_responses_ = function(dataSrc, qtpredicate=NULL, env=NULL, columns=c('person_id','item_id','item_score'))
{
  if(is.null(env))
    env = caller_env() 

  if(inherits(dataSrc,'data.frame'))
  {
    if(is.null(qtpredicate))
      return(dataSrc[,columns])    
    
    return(dataSrc[eval_tidy(qtpredicate, data=dataSrc, env=env), columns])
  }
  
  db_get_responses(dataSrc, qtpredicate=qtpredicate, env=env, columns=columns)
}


get_cte = function(db, columns)
{
  columns = setdiff(columns, 
    c(dbListFields(db,'dxResponses'), dbListFields(db,'dxScoring_rules')))
  
  cte = " dxResponses INNER JOIN dxScoring_rules USING(item_id, response)"
  
  if(length(intersect(dbListFields(db,'dxPersons'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxPersons USING(person_id)')

  if(length(intersect(dbListFields(db,'dxAdministrations'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxAdministrations USING(person_id, booklet_id)')
    
  if(length(intersect(dbListFields(db,'dxItems'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxItems USING(item_id)')
  
  if(length(intersect(dbListFields(db,'dxBooklets'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxBooklets USING(booklet_id)')
  
  if(length(intersect(dbListFields(db,'dxBooklet_design'),columns))>0) 
    cte = c(cte, 'INNER JOIN dxBooklet_design USING(booklet_id, item_id)')
  
  paste(cte, collapse=' ')
}

db_get_responses = function(db, qtpredicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) 
    env=caller_env() 
  
  columns = dbValid_colnames(columns)
  if(is.null(qtpredicate))
  {
    where = ''
    used_columns = columns
  } else
  {
    where = try(qtpredicate2where(qtpredicate, db, env), silent=TRUE)
    used_columns = union(columns, tolower(all.vars(qtpredicate)))                      
  }
  # decide which tables we need to join since joining costs time
  cte = get_cte(db, used_columns)
  
  if(!inherits(where, 'try-error'))
  { 
    respData = try(dbGetQuery(db, 
                     paste("SELECT", paste0(columns, collapse=','),
                           "FROM", cte, where)),
                   silent=TRUE)
    if(!inherits(respData,'try-error'))
      return(respData)
  }
  # to do: remove temp message
  message('sql translation did not work')
  print(qtpredicate)
  
  # translation to sql did not work  
  which_columns = intersect(c(columns, all.vars(qtpredicate)),
                            get_variables(db)$name)
          
  respData = dbGetQuery(db, 
               paste("SELECT", paste(which_columns, collapse=','),
                      "FROM", cte))

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
  if(is.null(qtpredicate))
  {
    where = ''
    used_columns = columns
  } else
  {
    where = try(qtpredicate2where(qtpredicate, db, env), silent=TRUE)
    # to do: non.nse??? 
    used_columns = union(columns, tolower(all.vars(qtpredicate)))                       
  }
  # decide which tables we need to join since joining costs time
  cte = get_cte(db, used_columns)
  
  if(!inherits(where, 'try-error'))
  { 
    # to do: postgres maes this a lower case column name, also for NFER its not an integer
    # use booklet_score everywhere, booklet_score is annoying
    respData = try(dbGetQuery(db, 
                     paste("SELECT", paste0(columns, collapse=','), ", SUM(item_score) AS booklet_score",
                             "FROM", cte, where,
                             "GROUP BY", paste0(columns, collapse=','),";")),        
                   silent=TRUE) 
    if(!inherits(respData,'try-error'))
      return(respData)
  }
  
  which_columns = intersect(c(columns, all.vars(qtpredicate),'item_score'),
                             get_variables(db)$name)
        
  respData = dbGetQuery(db, 
               paste("SELECT", paste(which_columns, collapse=','),
                        "FROM", cte))
  
  # this is a worse case scenario, it will take very long compared to using C, as a to do?
  return(
    respData[eval_tidy(qtpredicate, data=respData, env=env), columns] %>%
      group_by_at(columns) %>%
      summarise(booklet_score = sum(.data$item_score)) %>%
      ungroup())
}



# to do: better error message for sqlite if db connection is lost



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
            FROM dxAdministrations
              GROUP BY booklet_id)
            
         SELECT booklet_id, item_id, item_position, n_persons 
           FROM dxBooklet_design
             INNER JOIN bkl_count USING(booklet_id);"))
  
  if(is_bkl_safe(db, qtpredicate) )
  {
    uns = setdiff(c(dbListFields(db,'dxAdministrations'), dbListFields(db,'dxPersons')),
                  'booklet_id')
    
    where = try(qtpredicate2where(qtpredicate, db, env))
    
    if(!inherits(where,'try-error'))
    {
      ## to do: all.vars and sql
      if(length(intersect(all.vars(qtpredicate), uns)) == 0) 
      {
        res = try(
          dbGetQuery(db, paste(
            "WITH bkl_count AS(
              SELECT booklet_id, COUNT(*) AS n_persons 
                FROM dxAdministrations
                  GROUP BY booklet_id)
                
             SELECT booklet_id, item_id, item_position, n_persons 
               FROM dxBooklets
                 INNER JOIN bkl_count USING(booklet_id)
                   INNER JOIN dxBooklet_design USING(booklet_id)
                     INNER JOIN dxItems USING(item_id)",
              where,";")))
      } else
      {
        res = try(
          dbGetQuery(db, paste(
            "SELECT booklet_id, item_id, item_position, COUNT(*) AS n_persons 
               FROM dxBooklets
                 INNER JOIN dxBooklet_design USING(booklet_id)
                   INNER JOIN dxItems USING(item_id)
                    INNER JOIN dxAdministrations USING(booklet_id)
                      INNER JOIN dxPersons USING(person_id)",
                where,
            "GROUP BY booklet_id, item_id, item_position;")))
      }
      if(!inherits(res,'try-error'))
        return(res)
    }
  }
  
  # fail safe + case when !booklet_safe, will be slower  
  rsp = get_resp_data(db, qtpredicate, env=env, extra_columns = 'item_position')$x 
  if(nrow(rsp) == 0)
    return(rsp)
  
  rsp %>%
    group_by(.data$booklet_id, .data$item_id, .data$item_position) %>%
    summarise(n_persons = n()) %>%
    ungroup()

}

