df_identical = function(a, b)
{
  # check all values in dataframe equal, disregard column order
  
  if(!all(dim(a)==dim(b))) return(FALSE)
  if(!length(intersect(colnames(a),colnames(b))) == ncol(a)) return(FALSE)
  
  a = a %>% mutate_if(is.factor, as.character) 
  b = b %>% mutate_if(is.factor, as.character)
  
  for(col in colnames(a))
  {
    if(!all(a[,col]==b[,col])) return(FALSE)
  }
  return(TRUE)
}

dbRunScript <- function(db, fn)
{
  # run sql script included in the package
  # The R dbi api does not provide for execution of scripts.
  # A usual method to get around this is to split the input script on ;
  # however, there are numerous exceptions where this would not work (e.g. strings, triggers)
  # the kludge is to use a custom split string in our sql scripts, which is:
  # --#split#--

  fn = system.file("extdata", fn, package = "dexter", mustWork = TRUE)

  script = strsplit(paste0(readLines(fn, warn = FALSE), collapse='\n'),'--#split#--')[[1]]
  
  for (statement in script)
  {
    dbExecute(db,statement)
  }
}

# converts db from an old dexter version to the current dexter
convert_old_db = function(db)
{
  # assumes sqlite
  #dbRunScript(db,"dexter_sqlite.sql")
  dbTransaction(db,
  {
    if(dbExistsTable(db, "item_properties"))
    {
      iprops = setdiff(tolower(dbListFields(db,'item_properties')),'item')
      for(col in iprops)
      {
        dbExecute(db, paste0("ALTER TABLE dxItems ADD COLUMN ",col, sql_col_def('<empty>',is.default=TRUE),';'))
      }
      dbExecute(db, paste0('INSERT INTO dxItems(item_id,',paste(iprops,collapse=','),')
                              SELECT item,',paste(iprops,collapse=','),' FROM item_properties;'))
    } else 
    {
      dbExecute(db, 'INSERT INTO dxItems(item_id) SELECT item FROM Rules;')
    }
    
    dbExecute(db, 'INSERT INTO dxScoring_rules(item_id, response, item_score)
                    SELECT item, CAST(response AS TEXT), score FROM Rules;')
    if(dbExistsTable(db, "design"))
    {
      dbExecute(db, 'INSERT INTO dxBooklets(booklet_id) SELECT DISTINCT booklet_id FROM Design;')
      dbExecute(db, 'INSERT INTO dxTestparts(booklet_id, testpart_nbr) SELECT bookletName, 1 FROM Booklets;')
      dbExecute(db, 'INSERT INTO dxBooklet_design(booklet_id,	testpart_nbr,	item_id, item_position)
                        SELECT booklet_id, 1, item, position FROM Design;')
      dbExecute(db, 'INSERT INTO dxPersons(person_id) SELECT DISTINCT  dxpersonid FROM Persons;')
      dbExecute(db, 'INSERT INTO dxAdministrations(person_id, booklet_id) 
                        SELECT DISTINCT dxpersonid, booklet_id 
                          FROM Persons
                            INNER JOIN (SELECT DISTINCT booklet AS dxbookletid, booklet_id FROM Design) AS B1
                              USING(dxbookletid);')
      
      dbExecute(db, 'INSERT INTO dxResponses(person_id, booklet_id, item_id, response)
                        SELECT person, booklet_id, item, CAST(response AS TEXT)
                          FROM Responses
                            INNER JOIN (SELECT DISTINCT booklet, booklet_id FROM Design) AS B1
                                USING(booklet);')
      if('variable' %in% dbListFields(db,'persons'))
      {
        persons = dbGetQuery(db,'SELECT dxPersonID, variable, value FROM Persons;')
        persons %>% 
          group_by(.data$variable) %>%
          do({
            col = dbValid_colnames(.[1,1])
            dbExecute(db, paste0("ALTER TABLE dxPersons ADD COLUMN ",col, sql_col_def('<empty>',is.default=TRUE),';'))
            dbExecute(db, paste('UPDATE dxPersons SET',col,'=:val WHERE person_id=:person;'),
                      tibble(val = .$value, person = .$dxPersonID))
            data.frame()
          })
      }
    }
  })
}




                           
dbExists <- function(db, query, data)
{
  nrow(dbGetQuery(db, query, data)) > 0
}

dbUniquePersonIds <- function(db,n)
{
  if(is(db, 'SQLiteConnection'))
  {
    last = dbGetQuery(db,
                      "SELECT coalesce(MAX(CAST(substr(person_id,4) AS INTEGER)),0) AS n FROM dxPersons 
                        WHERE substr(person_id,1,3)='dxP' AND
                            upper(substr(person_id,4)) = lower(substr(person_id,4));")
  } 
  else if(is(db, 'PostgreSQLConnection')) 
  {
    last = dbGetQuery(db,"SELECT MAX(CAST(substring(person_id from 4) AS INTEGER)) 
                          FROM dxPersons WHERE person_id ~ '^dxP\\d+$';")
  }
 
  if (nrow(last)==0) { last = 0}  else {last = last[1,1]}
  
  return(paste0('dxP',c((last+1):(last+n))))
}

project_CreateTables <- function(db, covariates=NULL)
{
  # db is a handle to a database
  
  if(is(db, 'SQLiteConnection'))
  {
    dbRunScript(db,"dexter_sqlite.sql")
  } 
  else if(is(db, 'PostgreSQLConnection')) 
  {
    stop('Postgres is not supported yet')
    dbRunScript(db,"dexter_standard.sql")
    dbRunScript(db,"dexter_pgsql_triggers.sql")
  } 
  else {stop('unsupported database')}
  
  if (!is.null(covariates))
  {
    # do some cleaning to make sure these are acceptable column names
    names(covariates) = dbValid_colnames(names(covariates))
    covariates[['person_id']] = NULL
    for(col in names(covariates))
    {
      dbExecute(db, paste0("ALTER TABLE dxPersons ADD COLUMN ",col,sql_col_def(covariates[[col]],is.default=TRUE),';'))
    }
  }
}

sql_data_type = function(value)
{
  if(inherits(value,'Date')) return(' DATE ')
  else if(inherits(value,'factor')) return(' TEXT ')
  else if(inherits(value,'POSIXlt') | inherits(value,'POSIXt')) return(' DATETIME ')
  else if(typeof(value) == 'integer') return(' INTEGER ')
  else if(typeof(value) == 'double') return(' DOUBLE PRECISION ')
  else return(" TEXT ")
}

sql_col_def = function(value, is.default=FALSE, db=NULL)
{
  dt = sql_data_type(value)
  if(!is.default)
  {
    return(dt) 
  } else if(is.numeric(value))
  {
    return(paste(dt,'DEFAULT',value))
  } else
  {
    return(paste(dt,'DEFAULT',sql_quote(value,"'")))
  }
}


dbValid_colnames <- function(vec)
{
   gsub('^(?=\\d)','c',gsub('[^0-9a-z_]','_',tolower(vec)), perl=TRUE)
}

dbTransaction = function(db, expr, on_error = stop)
{
  if(is(db, 'SQLiteConnection')) dbExecute(db,'pragma foreign_keys=1;')
  dbBegin(db)
  tryCatch(expr, error=function(e){dbRollback(db);on_error(e);}, finally=NULL)
  dbCommit(db)
}

sql_IN = function(v)
{
  # returns IN(?[,?]*n) for vector v of length n
  paste0(' IN(',paste0(replicate(length(v),'?'),collapse=','),') ')
}


#' Selecting data
#' 
#' gather data from a dexter database 
#' 
#' @param dataSrc a dexter project database or data.frame
#' @param predicate an expression to select data on
#' @param columns the columns you wish to select
#' @param env optionally, an environment to evaluate the expression in
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
#' \dontrun{
#' # goal: fit the extended nominal response model using only persons 
#' # without any missing responses
#' library(dplyr)
#' 
#' # the following would not work since it will omit only the missing 
#' # responses, not the persons, which is not what we want in this case
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
#' }
get_responses <- function(dataSrc, predicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) env = caller_env() 
  get_responses_(dataSrc, eval(substitute(quote(predicate))), columns=columns, env=env)
}

#' Variables that can be used in expressions
#' 
#' Show the variables defined in your dexter project and their datatypes
#' 
#' @param db a dexter project database
#' @return a data.frame with name and type of the variables defined in your dexter project
#' @details 
#' A number of functions can take a dataSrc arguments and an optional expression. Expressions are a concise and
#' flexible way to filter data for the different psychometric functions in Dexter.
#' The variables that you can use in expressions consist of item properties and person covariates you defined
#' and a number of reserved variables that are automatically defined like \code{response} and \code{booklet_id}.
#' 
#' 
predicate_variables = function(db)
{
  cols = unique(
    lapply(c('dxItems','dxPersons','dxResponses','dxScoring_rules',
             'dxBooklets','dxBooklet_design'),
           function(tbl)
           {
             res = DBI::dbSendQuery(db,paste('SELECT * FROM',tbl,'WHERE 0=1;'))
             r = DBI::dbColumnInfo(res)
             DBI::dbClearResult(res)
             return(r)
           }) %>% 
      bind_rows()) %>% 
    arrange(.data$name)
  return(cols[cols$name != 'testpart_nbr',])
}


qtpredicate2where = function(qtpredicate, db, env)
{
  vars = unique(c(dbListFields(db,'dxItems'), dbListFields(db,'dxBooklets'), dbListFields(db,'dxPersons'), 
                  dbListFields(db,'dxBooklet_design'),dbListFields(db,'dxResponses')))
  sql = translate_sql_(list(partial_eval(qtpredicate,vars=vars,env=env)),con=db)
  
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
        #pff, strings. embedded quotes are doubled
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
  sql = gsub("IN *([^, \\(\\)'\"]+)","IN(\\1)", sql, perl=TRUE)
  sql = gsub("IN *('[^']*')","IN(\\1)", sql, perl=TRUE)
  sql = gsub('IN *("[^"]*")',"IN(\\1)", sql, perl=TRUE)
  
  if(length(sql) > 1) sql = paste0('(',sql,')',collapse=' AND ')
  
  return(paste(' WHERE ', sql))
}

get_responses_ <- function(dataSrc, qtpredicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) env=caller_env() 
  columns = dbValid_colnames(columns)
  if(inherits(dataSrc,'data.frame'))
  {
    if(is.null(qtpredicate))
    {
      respData = dataSrc[,columns]
    } else
    {
      respData = with(dataSrc, dataSrc[eval(qtpredicate), columns])
    }
  } else
  {
    db = dataSrc
    if(is.null(qtpredicate))
    {
      where = ''
      used_columns = columns
    } else
    {
      where = qtpredicate2where(qtpredicate, db, env)
      used_columns = union(columns, get.vars(qtpredicate))                      
    }
    # decide which tables we need to join since joining costs time
    used_columns = setdiff(used_columns, 
                           c(dbListFields(db,'dxResponses'), dbListFields(db,'dxScoring_rules')))
    cte = c()
    if(length(intersect(dbListFields(db,'dxPersons'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxPersons USING(person_id)')
    if(length(intersect(dbListFields(db,'dxItems'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxItems USING(item_id)')
    if(length(intersect(dbListFields(db,'dxBooklets'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxBooklets USING(booklet_id)')
    if(length(intersect(dbListFields(db,'dxBooklet_design'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxBooklet_design USING(booklet_id, item_id, testpart_nbr)')
    
    # can have unexpected behavior if person makes more than one testform
    respData = dbGetQuery(db, 
                          paste("SELECT", 
                                paste0(columns, collapse=','),
                                "FROM dxResponses INNER JOIN dxScoring_rules USING(item_id, response)",
                                paste0(cte,collapse=" "),
                                where,
                                'ORDER BY person_id, item_id;'))
  }
  return(respData)
}


is_bkl_safe = function(db, qtpredicate)
{
  # evaluates quoted expression used in _get_responses
  # to see if it is safe to trust the booklet from the database
  
  # we can be sure the booklets are not mutilated if
  # no item or respons level columns are used in the expression
  
  blacklist = unique(c( dbListFields(db,'dxItems'),
                        dbListFields(db,'dxScoring_rules'),
                        dbListFields(db,'dxBooklet_design'))) 
  
  blacklist = blacklist[blacklist!='booklet_id']

  return(length(intersect(get.vars(qtpredicate), blacklist)) == 0 )
}

get_sumscores <- function(dataSrc, qtpredicate=NULL,extra_columns=c(),env=NULL){
  # gets data and adds sumscores (and booklets if necessary)
  if(is.null(env)) env = caller_env()
  if(inherits(dataSrc,'data.frame')) { bkl_safe = FALSE
  } else if(is.null(qtpredicate)) { bkl_safe = TRUE 
  } else bkl_safe = is_bkl_safe(dataSrc, qtpredicate) 
  
  if(bkl_safe) { columns = c('person_id','item_id','item_score','booklet_id')
  } else columns = c('person_id','item_id','item_score')
  
  columns = union(columns, extra_columns)
  
  if(is.null(qtpredicate))
  {
    x = get_responses_(dataSrc, columns = columns) 
  } else
  {
    x = get_responses_(dataSrc, qtpredicate, columns, env = env) 
  }   
  # we now infer the design only if necessary
  if(!bkl_safe)
  {
    x = x %>%	arrange(.data$person_id, .data$item_id) 
    
    fitem = factor(x$item_id)
    x$iid = fmatch(x$item_id, fitem)
    # new dplyr requires ungroup
    x = x %>% 
      group_by(.data$person_id) %>% 
      mutate(sumScore = sum(.data$item_score), booklet_id=paste0(.data$iid, collapse=' ')) %>%
      ungroup()
    fbook = factor(x$booklet_id)
    x$booklet_id = fmatch(x$booklet_id, fbook)
  } else
  {
    x = x %>% 
      group_by(.data$person_id, .data$booklet_id) %>% 
      mutate(sumScore=sum(.data$item_score)) %>% 
      ungroup()
  }
  return(x)
}


get_design <- function(dataSrc, qtpredicate=NULL,env=NULL)
{
  # returns data.frame(booklet,item)
  # if you already have responsdata or sumscores, better use that to infer the design
  # to prevent drawing same data from the db twice
  
  if(is.null(env)) env = caller_env()
  if(inherits(dataSrc,'data.frame')) { bkl_safe = FALSE
  } else if(is.null(qtpredicate)) { bkl_safe = TRUE 
  } else bkl_safe = is_bkl_safe(dataSrc, qtpredicate) 
  
  if(!inherits(dataSrc,'data.frame') & bkl_safe) 
  {
    if(is.null(qtpredicate)){ where=''
    } else{ where = qtpredicate2where(qtpredicate, dataSrc, env)}
    
    return(dbGetQuery(dataSrc, 
                         paste('SELECT booklet_id, item_id, COUNT(*) AS n_persons
                                   FROM dxPersons
                                      INNER JOIN dxResponses USING(person_id)',
                                where,
                               'GROUP BY booklet_id, item_id;')))
  } else
  {
    return(get_sumscores(dataSrc, qtpredicate, env = env) %>% 
            group_by(.data$booklet_id, .data$item_id) %>%
            summarise(n_persons = n()) %>%
            ungroup() %>%
            select(.data$booklet_id, .data$item_id, .data$n_persons))
  }
}

get_design_from_sumscores = function(sumscores)
{
  sumscores %>%
    group_by(.data$booklet_id, .data$item_id) %>%
    summarise(n_persons = n()) %>%
    ungroup() %>%
    select(.data$booklet_id, .data$item_id, .data$n_persons)
}



