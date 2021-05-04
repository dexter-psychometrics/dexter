
is_db = function(db)
{
  inherits(db,'DBIConnection')
}

is_scored_db = function(db)
{
  !dbExists(db, 'SELECT 1 FROM dxScoring_rules WHERE CAST(item_score AS TEXT) <> response;')
}

dbRunScript = function(db, fn)
{
  fn = system.file("extdata", fn, package = "dexter", mustWork = TRUE)

  script = strsplit(paste0(readLines(fn, warn = FALSE), collapse='\n'),'--#split#--')[[1]]
  
  for (statement in script)
    dbExecute(db, statement)
}


dbExists = function(db, query, data=NULL){
  if(is.null(data))
    nrow(dbGetQuery(db, query)) > 0
  else
    nrow(dbGetQuery_param(db, query, data)) > 0
  
} 

dbCheck_reserved_colnames = function(nm)
{
  clash = intersect(tolower(nm), 
                    c('person_id','item_id','item_position',
                      'response','item_score','booklet_id'))
  
  if(length(clash) == 1)
  {
    stop(paste0("'", clash, "' is a reserved variable name in a dexter project"))
  } else if(length(clash) > 1)   
  {
    stop(paste(paste0("'",clash,"'", collapse=", "),
               'are reserved variable names in a dexter project'))
                
  }
}
  

dbUniquePersonIds = function(db,n)
{
  last = dbGetQuery(db,
                      "SELECT substr(person_id,4) AS n FROM dxpersons 
                        WHERE substr(person_id,1,3)='dx_' 
                          ORDER BY person_id DESC LIMIT 1;")

  if (NROW(last)==0) { last = 0L}  else {last = as.integer(last[1,1])}
  
  return(sprintf('dx_%07i',(1:n) + last))
}

project_CreateTables = function(db, person_properties=NULL)
{
  if(is(db, 'SQLiteConnection'))
  {
    dbRunScript(db,"dexter_sqlite.sql")
  } else
  {
    dbRunScript(db,"dexter_standard.sql")
  } 
  
  if (!is.null(person_properties))
  {
    # do some cleaning to make sure these are acceptable column names
    names(person_properties) = dbValid_colnames(names(person_properties))
    person_properties[['person_id']] = NULL
    for(col in names(person_properties))
    {
      dbExecute(db, paste0("ALTER TABLE dxpersons ADD COLUMN ",col,sql_col_def(person_properties[[col]],is.default=TRUE),';'))
    }
  }
}
# to do: depende on database, this is sqlite specific
sql_data_type = function(value)
{
  if(is.date(value))    return(' DATE ')
  if(is.factor(value))  return(' TEXT ')
  if(is.time(value))    return(' DATETIME ')
  if(is.integer(value)) return(' INTEGER ')
  if(is.numeric(value)) return(' DOUBLE PRECISION ')
  if(is.logical(value)) return(' INTEGER ')
  " TEXT "
}


sql_col_def = function(value, is.default=FALSE, db=NULL)
{
  dt = sql_data_type(value)
  
  if(!is.default || length(value)==0)
    return(dt) 
  
  if(is.date(value))
    return(paste0(dt," DEFAULT '",format(value, "%Y-%m-%d"),"'"))
  
  if(is.time(value))
    return(paste0(dt," DEFAULT '",format(value, "%Y-%m-%d %H:%M:%S"),"'"))

  if(is.numeric(value))
    return(paste(dt,'DEFAULT',if.else(is.na(value) || is.null(value),'NULL',value)))
  
  paste(dt,'DEFAULT',sql_quote(as.character(value),"'"))
  
}


dbValid_colnames = function(vec)
{
   gsub('^(?=\\d)','c',gsub('[^0-9a-z_]','_',tolower(vec)), perl=TRUE)
}

dbTransaction = function(db, expr, on_error = stop, on_error_rollback=TRUE)
{
  if(is(db, 'SQLiteConnection')) dbExecute(db,'pragma foreign_keys=1;')
  dbBegin(db)
  tryCatch(expr, error=function(e){if(on_error_rollback) dbRollback(db); on_error(e);}, finally=NULL)
  tryCatch(dbCommit(db), error=function(e){if(on_error_rollback) dbRollback(db); on_error(e);}, finally=NULL)
}


#to~do: for some reason data insertion in porstgres responses is extremely slow

# don't use literal strings containing : in these

dbGetQuery_param = function(db, statement, param)
{
  if(is(db, 'SQLiteConnection'))
    return(dbGetQuery(db,statement,param))
  
  param = as.list(param)
  
  if(is(db, 'PostgreSQLConnection') || is(db, 'PqConnection'))
  {
    vars = paste0(':',names(param))
    names(param) = NULL
    for(i in seq_along(vars))
    {
      statement = gsub(vars[i], paste0('$',i), statement, fixed=TRUE)
    }
  }  else if(is(db, 'RMySQL'))
  {
    vars = names(param)
    np = list()
    m = gregexpr('\\:\\w[\\w\\d_]*',statement)
    l = attr(m,'match.length')
    m = m[[1]]
    for(i in seq_along(m[[1]]))
    {
      if(m[i]>0)
      {
        var = substr(statement, m[i]+1, m[i]+l[i])
        if(var %in% vars)
        {
          statement = sub(var,'?',statement, fixed=TRUE)
          np[[length(np)+1]] = param[[var]]
        }
      }
    }
  }
  dbGetQuery(db,statement,param)
}
  
  

dbExecute_param = function(db, statement, param)
{
  if(is(db, 'SQLiteConnection'))
    return(dbExecute(db,statement,param))
  
  param = as.list(param)
  
  if(!endsWith(trimws(statement),';'))
    statement = paste0(statement,';')
  
  if(is(db, 'PostgreSQLConnection') || is(db, 'PqConnection'))
  {
    vars = paste0(':',names(param))
    names(param) = NULL
    for(i in seq_along(vars))
    {
      statement = gsub(paste0(vars[i],'(?=\\W)'), paste0('$',i), statement, perl=TRUE)
    }
  }  else if(is(db, 'RMySQL'))
  {
    vars = names(param)
    np = list()
    m = gregexpr('\\:\\w[\\w\\d_]*',statement)
    l = attr(m,'match.length')
    m = m[[1]]
    for(i in seq_along(m[[1]]))
    {
      if(m[i]>0)
      {
        var = substr(statement, m[i]+1, m[i]+l[i])
        if(var %in% vars)
        {
          statement = sub(var,'?',statement, fixed=TRUE)
          np[[length(np)+1]] = param[[var]]
        }
      }
    }
    param = np
  }
  
  dbExecute(db,statement,param)
}




