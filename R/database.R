
# change in rsqlite, now shows error if not all columns are used
# this is an ugly workaround, but it's a stupid change so hope it balances out
# not decided if should be used yet
dbGetQuery_ = function(db, stmt, df=NULL)
{
  if(is.null(df))
  {
    dbGetQuery(db, stmt)
  } else
  {
    m = gregexpr('(?<=[=,\\(]) *:\\w+', stmt, perl=TRUE)
    cols = gsub('^:','',trimws(regmatches(stmt, m)[[1]]), perl=TRUE)
    dbGetQuery(db, stmt, df[,cols])
  }
}

dbExecute_ = function(db, stmt, df=NULL)
{
  if(is.null(df))
  {
    dbExecute(db, stmt)
  } else
  {
    m = gregexpr('(?<=[=,\\(]) *:\\w+', stmt, perl=TRUE)
    cols = gsub('^:','',trimws(regmatches(stmt, m)[[1]]), perl=TRUE)
    dbExecute(db, stmt, df[,cols])
  }
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
  if(is(db, 'SQLiteConnection'))
  {
    dbRunScript(db,"dexter_sqlite.sql")
  } else if(is(db, 'PostgreSQLConnection')) 
  {
    stop('Postgres is not supported yet')
    #dbRunScript(db,"dexter_standard.sql")
    #dbRunScript(db,"dexter_pgsql_triggers.sql")
  } else {stop('unsupported database')}
  
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
  else if(inherits(value,'POSIXlt') || inherits(value,'POSIXt')) return(' DATETIME ')
  else if(typeof(value) == 'integer') return(' INTEGER ')
  else if(is.numeric(value)) return(' DOUBLE PRECISION ')
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
    return(paste(dt,'DEFAULT',ifelse(is.na(value) || is.null(value),'NULL',value)))
  } else
  {
    return(paste(dt,'DEFAULT',sql_quote(as.character(value),"'")))
  }
}


dbValid_colnames <- function(vec)
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

