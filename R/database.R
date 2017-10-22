

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

# converts db from a very old dexter version to the current dexter
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
    return(paste(dt,'DEFAULT',ifelse(is.na(value)|is.null(value),'NULL',value)))
  } else
  {
    return(paste(dt,'DEFAULT',sql_quote(value,"'")))
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

