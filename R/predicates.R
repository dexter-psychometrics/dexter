

# attempt to translate a quoted predicate and an environment to an SQL 'WHERE' statement
#
qtpredicate2where = function(qtpredicate, db, env)
{
  vars = unique(c(dbListFields(db,'dxItems'), dbListFields(db,'dxBooklets'), dbListFields(db,'dxPersons'), 
                  dbListFields(db,'dxBooklet_design'),dbListFields(db,'dxResponses')))
  sql = translate_sql_(list(partial_eval(qtpredicate,vars=vars,env=env)),con=db)

  if(getOption('dexter.debug', default=FALSE)) debug.log$send(sql, 'qtpredicate2where_uncorrected_sql')
  
  # translate_sql bug occurs for expressions like: booklet_id %in% as.character(1:4) 
  # <SQL> `booklet_id` IN CAST((1, 2, 3, 4) AS TEXT)
  # unfortunately regular expression functionality in R is a bit awkward 
  # so this is extremely kludgy to solve
  sql = lapply(sql, function(str)
  {
    split = strsplit(str,"CAST((", fixed=TRUE)[[1]]
    if(length(split)==1) return(str)
    split[1] = paste0(split[1],'(')
    split[2:length(split)] = vapply(split[2:length(split)], function(s)
    {
      mvec = regexpr('^.+(?=\\) AS )', s, perl=TRUE)
      vec = regmatches(s, mvec)
      remainder = substring(s,attr(mvec,'match.length')+6)
      mdtype = regexpr('^\\w+', remainder, perl=TRUE)
      dtype = regmatches(remainder, mdtype)
      remainder = substring(s, attr(mvec,'match.length') + 6 + attr(mdtype,'match.length')) 
      
      vec = trimws(vec)
      if(substr(vec,1,2) == "'")
      {
        #pff, strings. embedded quotes are already doubled
        vec = substr(vec, 2, nchar(vec)-1)			
        vec = paste0("CAST('", 
                     strsplit(vec,"(?<![^']')', ?'", perl=TRUE)[[1]],
                     "' AS ", dtype, ')', 
                     collapse=',') 
      } else
      {
        vec = paste0('CAST(', strsplit(vec, ',')[[1]], ' AS ', dtype, ')', collapse=',') 
      }
      paste0(vec, remainder)
    },"")
    paste0(split, collapse=' ')
  })
  

  # translate_sql_ has a bug with expressions like a %in% c(b) if b has length 1, solve here
  ## no longer necessary since dbplyr 1.1.0
  #sql = gsub("IN *([^, \\(\\)'\"]+)","IN(\\1)", sql, perl=TRUE)
  #sql = gsub("IN *('[^']*')","IN(\\1)", sql, perl=TRUE)
  #sql = gsub('IN *("[^"]*")',"IN(\\1)", sql, perl=TRUE)

  
  if(length(sql) > 1) sql = paste0('(',sql,')',collapse=' AND ')
  
  if(getOption('dexter.debug', default=FALSE)) debug.log$send(sql, 'qtpredicate2where_sql')
  
  return(paste(' WHERE ', sql))
}


# evaluates quoted expression as used in _get_responses
# to see if it is safe to trust the booklet_id column from the database
#  
# We can be sure the booklets are not mutilated if
# no item or respons level columns are used in the expression.
# This can err on the safe side but never on the fast but unsafe side.
is_bkl_safe = function(dataSrc, qtpredicate)
{
  if(inherits(dataSrc,'data.frame')) return(FALSE)
  if(is.null(qtpredicate)) return(TRUE)
 
  db = dataSrc
  
  blacklist = unique(c( dbListFields(db,'dxItems'),
                        dbListFields(db,'dxScoring_rules'),
                        dbListFields(db,'dxBooklet_design'))) 
  
  blacklist = blacklist[blacklist!='booklet_id']
  
  return(length(intersect(all.vars(qtpredicate), blacklist)) == 0 )
}




#' Variables that are defined in the project
#' 
#' Inspect the variables defined in your dexter project and their datatypes
#' 
#' @param db a dexter project database
#' @return a data.frame with name and type of the variables defined in your dexter project
#' @details 
#' The variables in Dexter consist of the item properties and person covariates you specified
#' and a number of reserved variables that are automatically defined like \code{response} and \code{booklet_id}.
#' 
#' Variables in Dexter are most useful when used in predicate expressions. A number of functions can 
#' take a dataSrc argument and an optional predicate. Predicates are 
#' a concise and flexible way to filter data for the different psychometric functions in Dexter.
#' 
#' The variables can also be used to retrieve data in \code{\link{get_responses}}
#' 
get_variables = function(db)
{
    lapply(c('dxItems','dxPersons','dxResponses','dxScoring_rules','dxBooklets','dxBooklet_design'),
           function(tbl)
           {
             res = DBI::dbSendQuery(db,paste('SELECT * FROM',tbl,'WHERE 0=1;'))
             r = DBI::dbColumnInfo(res)
             DBI::dbClearResult(res)
             return(r)
           }) %>% 
      bind_rows() %>% 
      distinct() %>%
      filter(.data$name != 'testpart_nbr') %>%
      arrange(.data$name)
}

