# dirty trick to avoid rewriting dplyr's do
# jennybc on https://github.com/STAT545-UBC/Discussion/issues/451
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

# prevent cran from complaining
booklet_id='_'

############################################
######      exported functions        ######
############################################


#' Start a new project
#'
#' Imports a complete set of scoring rules and starts a new project (data base)
#'
#'
#'
#' @param rules A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
#' The order is not important but spelling is. Any other columns will be ignored.
#' @param db A connection to an existing sqlite database or a string specifying a filename
#' for a new sqlite database to be created. If this name does not
#' contain a path, the file will be created in the work
#' directory. Any existing file with the same name will be overwritten.
#' @param covariates An optional list of person covariates. Names should correspond to covariates.
#' Values are treated as default values. The datatype will also be inferred from the values.
#' Known covariates will be imported (if supplied) in \code{add_booklet}. 
#' The datatypes should match those of your default values.
#' @return If the scoring rules pass a sanity check, a handle to the data base.
#' Otherwise, a data frame listing the problems found.
#' @details This package only works with closed items (e.g. likert, MC or possibly short answer): 
#' it does not score any open items.
#' The first step to creating a project is to import an exhaustive list of all items and
#' all admissible responses, along with the score that any of the latter will be given.
#' Responses may be integers or strings but they will always be treated as strings.
#' Scores must be integers, and the minimum score for an item must be 0.
#' When inputting data, all responses not specified in the rules can optionally be treated as
#' missing and ultimately scored 0, but it is good style to include the missing
#' responses in the list. NA values will be treated as the string 'NA'.
#'
#' @examples
#'\dontrun{
#' head(verbAggrRules)
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'    covariates=list(gender="<unknown>"))
#' }
#' 
start_new_project <- function(rules, db="dexter.db", covariates=NULL) {
  # for backward compatibility we rename if necessary
  names(rules)[names(rules)=='item'] = 'item_id'
  names(rules)[names(rules)=='score'] = 'item_score'
  
  rules = rules[, c("item_id", "response", "item_score")]
  rules$response = as.character(rules$response)
  rules$response[is.na(rules$response)] = 'NA'
  rules$item_id = as.character(rules$item_id)
  
  sanity = rules %>%
    group_by(.data$item_id) %>%
    summarise(less_than_two_scores = length(unique(.data$item_score))<2,
              duplicated_responses = any(duplicated(.data$response)),
              min_score_not_zero = min(.data$item_score)>0) %>%
    filter(.data$less_than_two_scores | .data$duplicated_responses | .data$min_score_not_zero)

  if (nrow(sanity)) {
    cat("There were problems with your scoring rules.\nCheck the output for possible reasons.\n")
    return(sanity)
  } else 
  {
    if (is(db,'character'))
    {
      if (file.exists(db)) file.remove(db)
      db = dbConnect(SQLite(), db)
    }
    dbTransaction(db,
    {
      project_CreateTables(db, covariates)
      dbExecute(db,'INSERT INTO dxItems(item_id) VALUES(:id);', 
                      tibble(id=unique(rules$item_id)))
      dbExecute(db,'INSERT INTO dxScoring_rules(item_id, response, item_score) 
                          VALUES(:item_id, :response, :item_score);', rules)
    })
    return(db)
  }
}



#' Open an existing project
#'
#' Opens a database created by function \code{start_new_project}
#'
#'
#' @param db_name The name of the data base to be opened.
#' @param convert_old automatically try to convert databases 
#' from older versions of Dexter that are no longer supported.
#' @return A handle to the data base.
#'
open_project <- function(db_name="dexter.db", convert_old=FALSE) {
  if (file.exists(db_name)) {
    db = dbConnect(SQLite(), db_name)
    dbExecute(db,'pragma foreign_keys=1;')
    if(!dbExistsTable(db,'dxItems'))
    {
      if(dbExistsTable(db,'rules'))
      {
        if(!convert_old){ 
          dbDisconnect(db)
          stop(paste('This appears to be a database from a previous version of Dexter',
                     'that is no longer supported. Use open_project with convert_old=TRUE',
                     'to attempt to convert your database to the newer version of Dexter.'))
        } else 
        { 
          cat('Converting database to the new Dexter version.')
          convert_old_db(db)
        }
      } else
      {
        dbDisconnect(db)
        stop('Sorry, this does not appear to be a Dexter database.')
      }
    }
  } else stop("There is no such file")
  return(db)
}


#' Add or modify scoring rules
#' 
#' Having to alter or add a scoring rule is occasionally necessary, e.g. in case of a key error. 
#' This function offers the possibility to do so and also allows you to add new items to your project
#' 
#' @param db handle to a Dexter project database
#' @param rules A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
#' The order is not important but spelling is. Any other columns will be ignored. See details
#' @return If the scoring rules pass a sanity check, a small summary of changes.
#' Otherwise, nothing.
#' @details 
#' The rules should contain all rules that you want to change or add. This means that in case of a key error
#' in a single multiple choice question, you typically have to change two rules.
#' @examples 
#'\dontrun{
#' # given that in your dexter project there is an mc item with id 'itm_01', 
#' # which currently has key 'A' but you want to change it to 'C'.
#' 
#' new_rules = data.frame(item_id='itm_01', response=c('A','C'), item_score=c(0,1))
#' touch_rules(db, new_rules)
#' }
#' 
touch_rules = function(db, rules)
{
  # for backward compatibility we rename if necessary
  names(rules)[names(rules)=='item'] = 'item_id'
  names(rules)[names(rules)=='score'] = 'item_score'
  
  rules = rules[, c("item_id", "response", "item_score")]
  rules$response = as.character(rules$response)
  rules$response[is.na(rules$response)] = 'NA'
  # the following line gets rid of dplyr factor warnings
  rules$item_id = as.character(rules$item_id)
  
  existing_rules = dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxScoring_rules;')
  existing_opts = existing_rules %>% select(-.data$item_score)
  
  new_rules = rules %>% anti_join(existing_opts, by=c('item_id','response'))
  amended_rules = rules %>% inner_join(existing_opts, by=c('item_id','response'))
  
  # to judge the validity of the new rules, we have to look at them in combination
  # withj the rules in the db that will not be changed
  sanity = new_rules %>% 
    dplyr::union(amended_rules)  %>% 
    dplyr::union(existing_rules %>% 
                   inner_join(tibble(item_id=amended_rules$item_id), by='item_id') %>%
                   anti_join(amended_rules, by=c('item_id','response'))
                 ) %>%
    group_by(.data$item_id) %>%
    summarise(less_than_two_scores = length(unique(.data$item_score))<2,
              duplicated_responses = any(duplicated(.data$response)),
              min_score_not_zero = min(.data$item_score)>0) %>%
    filter(.data$less_than_two_scores | .data$duplicated_responses | .data$min_score_not_zero)
  
  
  if (nrow(sanity)) {
    cat("There were problems with your scoring rules.\nCheck the output for possible reasons.\n")
    return(sanity)
  }    
  
  dbTransaction(db,
  {
    if(nrow(new_rules)>0)
    {
      new_items = setdiff(new_rules$item_id, dbGetQuery('SELECT item_id FROM dxItems;')$item_id)
      if(length(new_items)>0) dbExecute(db,'INSERT INTO dxItems(item_id) VALUES(:id);', tibble(id=new_items))
      dbExecute(db,'INSERT INTO dxScoring_rules(item_id, response, item_score) 
							                            VALUES(:item_id, :response, :item_score);', new_rules)
    }
    if(nrow(amended_rules)>0) 
    {
      dbExecute(db,'UPDATE dxScoring_rules SET item_score=:item_score 
                      WHERE item_id=:item_id AND response=:response;', 
                amended_rules)
    }
  })
  cat(paste0('rules_changed: ', nrow(amended_rules), '\nrules_added: ', nrow(new_rules)))
}

#' Show scoring rules
#' 
#' Show the scoring rules currently present in the dexter project db
#' 
#' @param db handle to a Dexter project database
#' @return data.frame of scoring rules
#' 
show_rules = function(db)
{
  dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxScoring_rules ORDER BY item_id, response;')
}


#' Add a booklet to a project
#'
#' Adds item response data for a test form (a.k.a. booklet)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param x A data frame containing the responses and, possibly, some additional
#' person characteristics. See details.
#' @param booklet_id A (short) string identifying the test form (booklet)
#' @param auto_add_unknown_rules  If FALSE, an error will be generated if 
#' some of the responses do not appear in the scoring rules. Default is TRUE.
#' @return A list of: \item{items}{The names of the columns in \code{x} that were
#' treated as items}
#' \item{covariates}{The names of the columns in \code{x} that were
#' treated as person covariates}
#' \item{not_listed}{A data frame of all responses that will be treated as missing}
#' @details It is common practice to keep data in rectangular tables: data frames
#' or foreign software like Excel, SPSS, etc. This function is provided to input
#' data in that form, one booklet at a time. The starting point is a data frame,
#' and getting the data frame into R is left to the user. We have found package
#' \code{readxl} to be very good at reading Excel sheets, and \code{haven} quite
#' efficient with SPSS files.
#'
#' If the dataframe \code{x} contains a variable named \code{person_id} this variable 
#' will be used to identify unique persons. It is assumed that a single person will only 
#' make a single booklet once, otherwise an error will be generated. 
#' 
#' If a person ID is not supplied, dexter will generate unique person ID's for each row of data.  
#'
#' Any variable whose name has an exact match in the scoring rules input with
#' function \code{start_new_project} will be treated as an item; any variable whose name has an 
#' exact match in the covariates will be treated as covariate. If a name matches both
#' a covariate and an item, the item takes precedence. Variables other than items, covariates 
#' and person_id will be ignored.
#' 
#' If \code{auto_add_unknown_rules=TRUE}, any responses to an item that do not have an 
#' exact match in the scoring rules will be treated 
#' as missing  and ultimately given the lowest score of 0. 
#' To score missing data differently, 
#' or simply abide to good style, the user can include explicit entries for missing value
#' indicators in the scoring rules.
#' 
#' @examples 
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'   covariates=list(gender="<unknown>"))
#' head(verbAggrData)
#' add_booklet(db, verbAggrData, "agg")
#' }
#' 
add_booklet <- function(db, x, booklet_id, auto_add_unknown_rules = TRUE) {
  x = x %>% mutate_if(is.factor, as.character) 
  
  covariates = intersect(dbListFields(db,'dxPersons'), tolower(names(x)))
  covariates = covariates[covariates != 'person_id']
  
  design = tibble(booklet_id = booklet_id, item_id = names(x), col_order=c(1:ncol(x))) %>%
    inner_join(dbGetQuery(db, "SELECT item_id FROM dxItems;"), by='item_id') %>%
    arrange(.data$col_order) %>% 
    select(-.data$col_order)
  
  design$item_position = c(1:nrow(design))

  dbTransaction(db,{
    if (dbExists(db,'SELECT 1 FROM dxBooklets WHERE booklet_id=:b;',tibble(b=booklet_id)) ) 
    {
      if(!df_identical(dbGetQuery(db,'SELECT item_id
                                      FROM dxBooklet_design
                                      WHERE booklet_id=:b
                                     ORDER BY testpart_nbr,item_position;', tibble(b=booklet_id)),
                    tibble(item_id = design$item_id)))
      {
        stop("There is already a booklet with this ID which has different items or a different item order")
      }
    } else
    {  
      dbExecute(db,'INSERT INTO dxBooklets(booklet_id) VALUES(:b);',tibble(b=booklet_id))
      dbExecute(db,'INSERT INTO dxTestparts(booklet_id, testpart_nbr) VALUES(:b,:tp);',tibble(b=booklet_id,tp=1))
      dbExecute(db,'INSERT INTO dxBooklet_design(booklet_id,testpart_nbr, item_id, item_position) 
                          VALUES(:booklet_id,1,:item_id,:item_position);', design)
    }
                  
    x$booklet_id = booklet_id
    if(!'person_id' %in% names(x))
    {
      x$person_id = dbUniquePersonIds(db,nrow(x))
      new_people = x$person_id
    } else
    {
      known_people = dbGetQuery(db,'SELECT person_id FROM dxPersons;')$person_id
      new_people = setdiff(x$person_id,known_people)
    }
                  
    dbExecute(db,'INSERT INTO dxPersons(person_id) VALUES(:person_id);', tibble(person_id=new_people))
    dbExecute(db,'INSERT INTO dxAdministrations(person_id,booklet_id) VALUES(:person_id,:b);', tibble(person_id=new_people,b=booklet_id))
          
    responses = x[,c(design$item_id, "booklet_id", "person_id")] %>%
      gather_(key_col='item_id', value_col='response', gather_cols=design$item_id, na.rm=FALSE)
            
    responses$response = as.character(responses$response)
    responses$response[is.na(responses$response)] = 'NA'
                  
    if(auto_add_unknown_rules)
    {
      existing_rules = dbGetQuery(db, "SELECT item_id, response FROM dxScoring_rules;")
      rules = responses[,c('item_id', 'response')]
      new_rules = rules[!duplicated(rbind(existing_rules, rules))[-seq_len(nrow(existing_rules))], ]
      if (nrow(new_rules)>0) dbExecute(db,'INSERT INTO dxScoring_rules(item_id,response,item_score) VALUES(:item_id,:response,0);',new_rules)
    } else
    {
      new_rules=NA
    }
    dbExecute(db,'INSERT INTO dxResponses(booklet_id,person_id,item_id,response) 
                                VALUES(:booklet_id,:person_id,:item_id,:response);', responses)
            
    # make this report before we mutilate the colnames  
    columns_ignored = setdiff(names(x), c(design$item_id,'person_id','item_id','booklet_id') )
    columns_ignored = columns_ignored[!tolower(columns_ignored) %in% covariates]
                  
    # add covariates
    if(length(covariates)>0)
    {
      names(x) = tolower(names(x))
      dbExecute(db,paste0('UPDATE dxPersons SET ',paste0(covariates,'=:',covariates,collapse=','),' WHERE person_id=:person_id;'),
                               x[,c(covariates,'person_id')])
    }
  }) 
  
  return(
    list(
      items = design$item_id,
      covariates = covariates,
      columns_ignored = columns_ignored,
      auto_add_unknown_rules=auto_add_unknown_rules,
      zero_rules_added = new_rules
    )
  )
 }



#' Add item properties to a project
#'
#' Adds item properties to an existing database
#'
#'
#' @param db A handle to the database, e.g. the output of \code{start_new_project}
#' or \code{open_project}
#' @param item_properties A data frame containing the item properties. See details.
#' @param overwrite Whether overwrite is permitted (default=FALSE)
#' @return A list of: \item{unknown_items}{Item IDs for any items that were provided
#' in the data frame but could not be found in the data base}
#' \item{items_unaccounted_for}{Item IDs for any items that exist in the data base
#' but were not given properties in the data frame}
#' @details When entering response data in the form of a rectangular person x item
#' table, it is easy to provide person properties but practically impossible
#' to provide item properties. This function provides a possibility to do so.
#' The order of the rows and columns in the data frame is not important but
#' (i) there must be a column called exactly \code{item} or \code{item_id} containing the item IDs
#' exactly as entered before, and (ii) all items in the data frame must be known
#' to the data base and all items in the data base must be given properties --
#' otherwise, there will be a warning message, and nothing else will be done.
#' If all is well, the data frame will be added to the project database, and any variables in 
#' it may be used in analyses involving item properties.
#' 
#' @seealso \code{\link{fit_domains}}, \code{\link{profile_plot}} for
#'  possible uses of item_properties
#'
#' @examples 
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'   covariates=list(gender="<unknown>"))
#' head(verbAggrProperties)
#' add_item_properties(db, verbAggrProperties)
#' show_item_properties(db)
#' }
#' 
#'
add_item_properties <- function(db, item_properties, overwrite=FALSE) {
  item_properties = item_properties %>%
    mutate_if(is.factor, as.character) 
  
  names(item_properties)[names(item_properties)=='item'] = 'item_id'
  names(item_properties) = dbValid_colnames(names(item_properties))
  
  if(!'item_id'%in%names(item_properties))
  {
    stop("there was no column provided with name 'item_id'")
  }
  existing_item_properties = dbListFields(db, 'dxItems') # for convenience we include the item_id as a property
  if(!overwrite & !setequal(intersect(names(item_properties), existing_item_properties), 'item_id'))
  {
    stop('Some of the listed item properties already exist, specify overwrite=TRUE to overwrite')
  }

  if(!setequal(dbGetQuery(db,'SELECT item_id FROM dxItems;')$item_id, item_properties$item_id))
  {
    stop('properties not specified for all items.')
  }
  dbTransaction(db, 
  {
    for(prop_name in setdiff(names(item_properties), existing_item_properties))
    {
      dbExecute(db, paste0("ALTER TABLE dxItems ADD COLUMN ",prop_name, sql_data_type(item_properties[,prop_name]),";"))
    }
    pnames = names(item_properties)[names(item_properties)!='item_id']
    
    dbExecute(db,paste0('UPDATE dxItems SET ',paste0(pnames,'=:',pnames,collapse=', '),' WHERE item_id=:item_id;'),
               item_properties)
  })
  invisible(NULL)
}




#' List booklets in a project
#'
#' Show a list of the test forms (booklets) that have been entered in the db
#' so far
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame showing the booklet_id, the number of persons and the number of items.
#'
show_booklets <- function(db) {
  dbGetQuery(db,'SELECT booklet_id, n_items, n_persons FROM dxBooklet_stats ORDER BY booklet_id; ')
}


#' List item properties
#'
#' Show a list of the item properties defined in the project (if any)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data.frame of item properties, values and counts
#'
show_item_properties <- function(db) {
  dbGetQuery(db, 'SELECT * FROM dxItems;') %>% 
    select(-.data$item_id) %>% 
    map_df(function(x){tibble(value=as.character(x)) %>% group_by(.data$value) %>% mutate(N=n())}, .id='item_property')
}



#' List person properties
#'
#' Show all person properties defined in the project (if any)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data.frame of person properties, values and counts
#'
show_person_properties <- function(db) {
  dbGetQuery(db, 'SELECT * FROM dxPersons;') %>% 
    select(-.data$person_id) %>% 
    map_df(function(x){tibble(value=as.character(x)) %>% group_by(.data$value) %>% mutate(N=n())}, .id='person_property')
}



#' List items in a project
#'
#' Show all items that have been entered in the db
#' so far together with the item properties
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with items and item properties
#'
show_items <- function(db){
  dbGetQuery(db,'SELECT * FROM dxItems ORDER BY item_id;')
}

#' Show the test design
#'
#' Show a list of all items that have been entered in the db
#' so far by booklet and position in the booklet
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame (in string format) showing the booklet design: rows are items,
#' columns are booklets, and the numbers in the cells show the position of the item in the
#' booklet (blank if the booklet does not include the item)
#'
show_design <- function(db){
  design = dbGetQuery(db,'SELECT booklet_id, item_id, item_position FROM dxBooklet_design;')
  bar = format(spread_(design, key_col='booklet_id', value_col='item_position'))
  bar[bar=="NA"] = ""
  bar
}

#' Show persons in a project
#'
#' Show all persons that have been entered in the db
#' so far together with covariates
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with persons and covariates
#'
show_persons <- function(db){
  dbGetQuery(db,'SELECT * FROM dxPersons ORDER BY person_id;')
}

#' Provide test scores
#'
#' Supplies the weighted sum of item scores for each person selected.
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @return A data frame showing person_id' s and their test scores
#' 
get_testscores<-function(dataSrc, predicate=NULL) {
  qtpredicate=eval(substitute(quote(predicate)))
  x = get_sumscores(dataSrc, qtpredicate, env=caller_env())
  if(nrow(x) == 0) stop('empty selection')
  x %>% group_by(.data$person_id) %>% 
    slice(1) %>% 
    ungroup() %>%
    select(.data$person_id, .data$booklet_id, test_score=.data$sumScore) 
}


#' Design as network
#'
#' Export the test design as an incidence matrix 
#' and a weight matrix
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param weights Weight the edges between booklets by the number of common 
#' \code{"items"} or \code{"responses"} (default is items). 
#' @return A list of two data frames: \item{im}{incidence matrix}
#' \item{wm}{weights matrix}
#' @details 
#' The output of this function can be passed to packages for network analysis
#' such as \code{igraph} or \code{qgraph}. We prefer to not load these 
#' packages automatically as they are fairly large and rely on a number 
#' of dependencies. 
#' @examples
#' \dontrun{
#' dsgn = design_as_network(db)
#' # Check if design is connected
#' design_is_connected(dsgn)
#' }
#'
design_as_network <- function(dataSrc, predicate = NULL, weights=c("items","responses")){
  w = match.arg(weights)
  qtpredicate = eval(substitute(quote(predicate)))
  
  design = get_design(dataSrc, qtpredicate, env=caller_env())
  
  if(length(unique(design$booklet_id)) < 2) stop("This makes sense only if you have at least two booklets")
  
  im = as.matrix(table(design$item_id, design$booklet_id))
  wm = crossprod(im, im)
  diag(wm) = 0
  if (w=="responses") {
    b = design %>% group_by(.data$booklet_id) %>% slice(1)
    ww = outer(b$n_persons, b$n_persons, "+")
    wm = wm*ww    
  }
  list(im=im, wm=wm)
}  

#' Test if design is connected
#'
#' Use the output from design_as_network to check if your design is connected.
#'
#'
#' @param design Output from design_as_network
#' @return TRUE or FALSE
#' @examples
#' \dontrun{
#' # as an example, turn off some your booklets and see if you are
#' # still left with a connected design
#' dsgn = design_as_network(db, !(booklet_id %in% c('b1','b3','b4')))
#' design_is_connected(dsgn)
#' }
#'
design_is_connected = function(design)
{
  # implementation of depth first search
  d = design$wm
  visited = rep(FALSE, ncol(d))
  rownames(d) = c(1:nrow(d))
  colnames(d) = c(1:nrow(d))
  dfs = function(start)
  {
    start = as.integer(start)
    if(visited[start]) return(0)
    visited[start] <<- TRUE
    vapply(rownames(d)[d[,start]>0], dfs, 0)
    0
  }
  dfs(1)
  return(all(visited))
}


#' Estimate the Rasch and the Interaction model
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @return An object of class \code{rim} holding results
#' for the Rasch model and the interaction model.
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function fits the
#' Rasch model and the interaction model on a complete rectangular array of
#' responses, with comparison between the two models playing an important role.
#' The rectangular array is typically one booklet but can also consist of 
#' the intersection (common items) of two or more booklets. If the intersection is empty 
#' (no common items for all persons), the function will exit with an error message.
#'
#' @seealso \code{\link{plot.rim}}, \code{\link{fit_domains}} 
#' 
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' 
#' m = fit_inter(db, booklet_id=='agg')
#' plot(m, "S1DoScold", show.observed=TRUE)
#' }
fit_inter <- function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  respData = get_responses_(dataSrc,qtpredicate, env=caller_env()) 
  
  if(is.null(qtpredicate)){ grpName = 'All'} else{ grpName = as.character(qtpredicate)}

  # make the sufficient stats
  # for this function a subset is required which all students make
  # this is actually easier to do with wide format  
  if (nrow(respData)<1) stop("No responses to analyse")
  
  respData = 
    respData %>%
    spread_(key_col='item_id', value_col='item_score') %>%
    select(which(c(TRUE, !is.na(colSums(.[,-1]))))) 
  
  if(ncol(respData)==1) stop (paste('The instersection of responses in your data is empty.',
                                    'The interaction model cannot be computed concurrently for a whole design of test forms.',
                                    'Type ?fit_inter for more informtation'))
    
  # derive items to add zero scores, reportedly not necessary in future
  score0 = tibble(item_id = names(respData)[2:ncol(respData)], item_score=0)
    
  itcol = setdiff((names(respData)), 'person_id')
  respData = respData %>%
    gather_(key_col='item_id', value_col='item_score', itcol)

  respData = respData %>% 
           group_by(.data$person_id) %>% 
           mutate(sumScore=sum(.data$item_score))
  
  ssScoreLev = respData %>% 
           group_by(.data$item_id, .data$item_score) %>% 
           summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>% 
           full_join(score0, by = c("item_id","item_score")) %>%
           arrange(.data$item_id, .data$item_score)
  ssScoreLev[is.na(ssScoreLev)] = 0
  # here we have to add ungroup because the new dplyr syntax otherwise won't work
  ssItemLev = ssScoreLev %>% ungroup() %>%
           group_by(.data$item_id) %>%
           summarise(nCat = n(), N = sum(.data$sufI), sufC = sum(.data$sufC)) %>%
           mutate(first = cumsum(.data$nCat) - .data$nCat + 1, last = cumsum(.data$nCat))  %>%
           arrange(.data$item_id)
  plt = respData %>%
           group_by(.data$item_id,.data$sumScore) %>% 
           summarise(meanScore = mean(.data$item_score), N = n())
  maxTotScore = sum(tapply(ssScoreLev$item_score, ssScoreLev$item_id, max))
  totalLev = plt[plt$item_id==plt$item_id[1],] %>% 
           right_join(tibble(sumScore=0:maxTotScore), by="sumScore")
  totalLev$N[is.na(totalLev$N)] = 0
  ss = list(group=grpName, il=ssItemLev, sl=ssScoreLev, tl=totalLev, plt=plt)
  result = try(EstIM(ss))
  if (inherits(result, "try-error")) result=NULL
    # add the regressions, convenient for plotting
  if (!is.null(result)) 
  {
    C = rep(1:nrow(ss$il), ss$il$nCat)
    ctrRM=ittotmat(result$bRM, result$cRM[C], ss$sl$item_score, ss$il$first, ss$il$last)
    ctrIM=ittotmat(result$bIM, result$cIM[C], ss$sl$item_score, ss$il$first, ss$il$last)
    mm = sweep(model.matrix(~0+ss$sl$item_id), 1, ss$sl$item_score, '*')
    itrRM = as.data.frame(crossprod(mm, ctrRM))
    itrIM = as.data.frame(crossprod(mm, ctrIM))
    row.names(itrRM) = row.names(itrIM) = ss$il$item_id
    regs = list(ctrRM, ctrIM, itrRM, itrIM)
    names(regs) = c('ctrRM','ctrIM','itrRM','itrIM')
  } else {  
    regs = NULL
  }
    outpt = list(est=result, ss=ss, regs=regs)
    class(outpt) = append("rim",class(outpt))
    outpt
}



###################################################
#' Distractor plot
#'
#' Produce a diagnostic distractor plot for an item
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param item The ID of the item to plot. A separate plot will be produced
#' for each booklet that contains the item, or an error message if the item ID
#' is not known. Each plot contains a non-parametric regression of each possible
#' response on the total score.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' 
distractor_plot <- function(dataSrc, item, predicate = NULL, nc=1, nr=1){  
  qtpredicate = eval(substitute(quote(predicate)))
  if(!inherits(dataSrc,'data.frame'))
  {
    if(is.null(qtpredicate))
    {
      #attempt to speed up a little
      booklets = dbGetQuery(dataSrc,'SELECT booklet_id FROM dxBooklet_design WHERE item_id=:item;',tibble(item=item))$booklet_id
      qtpredicate = quote(booklet_id %in% booklets)
      x = get_sumscores(dataSrc, qtpredicate, extra_columns=c('response','item_position'), env=environment())
    } else 
    {
      x = get_sumscores(dataSrc, qtpredicate, extra_columns=c('response','item_position'), env=caller_env())
    }
  } else 
  {
    x = get_sumscores(dataSrc, qtpredicate, extra_columns=intersect(c('response','item_position'),names(dataSrc)), env=caller_env())
  }
  x = x[x$item_id==item,]
  if(!'item_position' %in% names(x)) x$item_position = 0
  booklets = x %>%  
    select(.data$booklet_id, .data$item_position) %>% 
    group_by(.data$booklet_id) %>% 
    slice(1) %>% 
    ungroup() 
  if (nrow(x) < 1) stop(paste("Item", item, "not found in dataSrc."))
  names(x)[names(x)=='sumScore'] = 'booklet_score'

  foo = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$response, .data$item_score, .data$booklet_score) %>%
    summarise(n = n()) %>% ungroup()
  
  iSt = x %>%
    group_by(.data$booklet_id, .data$item_id) %>% 
    summarise(meanScore = mean(.data$item_score), 
              rit = cor(.data$item_score, .data$booklet_score), 
              rir = cor(.data$item_score, .data$booklet_score - .data$item_score)) %>% 
    inner_join(booklets, by = c("booklet_id"))
  
  mxSc = max(x$item_score)
  iSt$pvalue = iSt$meanScore / mxSc
  
  npic = nrow(booklets)
  ly = my_layout(npic, nr, nc)
  graphics::layout(matrix(1:(ly$nr * ly$nc), byrow = TRUE, 
                          ncol = ly$nc))
  
  foo$response=factor(foo$response)
  foo$response=addNA(foo$response, ifany=TRUE)
  labelz = levels(foo$response)
  
  lapply(iSt$booklet_id, function(x){
    st = as.list(iSt[iSt$booklet_id==x,])
    tit = sprintf("%s: position %d in %s", st$item_id, st$item_position, st$booklet_id)
    subtit = sprintf("Pval: %.3f, Rit: %.3f, Rir: %.3f", st$pvalue, st$rit, st$rir)
    y = foo[foo$booklet_id == x,]
    
    graphics::plot(c(0, max(y$booklet_score)), c(0, 1), type = "n", 
                   main = tit, sub = subtit, xlab = "Sum score", ylab = "Proportion", 
                   cex.sub = 0.8)
    
    bar = y %>% 
      group_by(.data$booklet_score) %>% 
      summarise(n = sum(.data$n))
    
    dAll = density(bar$booklet_score, n = 51, weights = bar$n/sum(bar$n))
    nnn = sum(bar$n)
    lgnd = y %>% group_by(.data$response)  %>% do({
      dxi = density(.$booklet_score, n = 51, weights = .$n/sum(.$n), 
                    bw = dAll$bw, from = min(dAll$x), to = max(dAll$x))
      yy = dxi$y/dAll$y * sum(.$n)/nnn
      k = match(.$response[1], labelz) + 1
      graphics::lines(dAll$x, yy, co = k, lw = 2)
      tibble(col = k, resp = paste0(.$response[1]," (", .$item_score[1], ")"))
    })
     graphics::legend("right", legend = as.character(lgnd$resp), 
                       lty = 1, col = lgnd$col, cex = 0.8, box.lty = 0)
      graphics::box()
   })
  return(NULL)
}



#' Estimate the Rasch and the Interaction model per domain
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param item_property The item property defining the 
#' domains (subtests)
#' @return An object of class \code{imp} holding results
#' for the Rasch model and the interaction model.
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function fits the
#' Rasch model and the interaction model on a complete rectangular array of
#' responses, with comparison between the two models playing an important role.
#' The rectangular array is typically one booklet but can also consist of 
#' the intersection (common items) of two or more booklets. If the intersection is empty 
#' (no common items for all persons), the function will exit with an error message.
#'
#' @seealso \code{\link{plot.rim}}, \code{\link{fit_inter}}, \code{\link{add_item_properties}} 
#' 
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' mSit = fit_domains(db, item_property= "situation")
#' plot(mSit)
#' }
#' 
fit_domains = function(dataSrc, item_property, predicate = NULL)
{
  columns = c('person_id','item_id','item_score', item_property)
  qtpredicate = eval(substitute(quote(predicate))) #magic!
  respData = get_responses_(dataSrc,qtpredicate, columns=columns,env=caller_env()) 
  if(nrow(respData) == 0) stop('no data to analyse')
  
  respData = respData %>% 
      group_by(.data$person_id, .data[[!!item_property]]) %>% 
      summarise(domain_score=sum(.data$item_score)) %>%
      select(.data$person_id, .data$domain_score, .data[[!!item_property]]) %>%
      ungroup()
  
  names(respData)[names(respData)==item_property] = 'item_id'
  names(respData)[names(respData)=='domain_score'] = 'item_score'
  
  fit_inter(respData)
}


#' Profile plot
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param item_property The item property defining the domains
#' @param covariate Person covariate name for subsetting
#' @param model "IM" (default) or "RM" where "IM" is the interaction model and 
#' "RM" the Rasch model. The interaction model is the default as it fits 
#' the data better at least as good as the Rasch model.
#' @param x Which value of the item_property to draw on the x axis, if NULL, one is chosen automatically
#' @param ... further arguments to plot, many have useful defaults
#' @return Nothing interesting
#' @details 
#' The profile plots can be used to investigate whether two (or more) groups of respondents 
#' attain the same test score in the same way. The user must provide a  
#' (meaningfull) classification of the items in two non-overlapping subsets such that 
#' the test score is the sum of the scores on the subsets. 
#' The plot shows the probabilities to obtain 
#' any combinations of subset scores with thin gray lines indicating the combinations 
#' that give the same test score. The thick lines connect the most likely 
#' combination for each test score in each group.
#' When applied to educational test data, the plots can be used to detect differences in the 
#' relative difficulty of (sets of) items for respondents that belong to different 
#' groups and are matched on the test score. This provides a content-driven way to 
#' investigate differential item functioning. 
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'    covariates=list(gender="<unknown>"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' profile_plot(db, item_property='mode', covariate='gender')
#' }
#' 
profile_plot <- function(dataSrc, item_property, covariate, predicate = NULL, model = "IM", x = NULL, ...) 
{
  if (model != "IM") model="RM"
  user.args = list(...)
  if(!inherits(dataSrc,'data.frame'))
  {
    item_property = tolower(item_property)
    covariate = tolower(covariate)
  }
  
  columns = c('person_id','item_id','item_score',item_property, covariate)
  
  qtpredicate=eval(substitute(quote(predicate)))
  respData = get_responses_(dataSrc, qtpredicate, columns = columns, env = caller_env()) 
  if(nrow(respData) == 0) stop('no data to analyse')
  
  prop_design = unique(respData[,c(item_property,'item_id')])
  props = unique(prop_design[[item_property]])
  if(length(props) != 2)
    stop('this function needs an item_property with 2 unique values in your data')
  if(!is.null(x)) if(x %in% props) props = c(x, props[props!=x])

  # split into list
  respData = by(respData, respData[[covariate]], identity)

  common_scores = Reduce(dplyr::intersect, lapply(respData, function(x) x[,c('item_id','item_score')]))
  problems = Reduce(dplyr::union,
                    lapply(respData,
                           function(data) dplyr::setdiff(data[,c('item_id','item_score')],common_scores)))
 
  if(nrow(problems) > 0)
  {
    problems = tibble(item_id = unique(problems$item_id))
    warning(paste('the following items do not have the same score categories over all covariates and',
                  'have been removed from the analysis:', paste0(problems$item_id,collapse=', ')))
    respData = lapply(respData, function(x) x %>% anti_join(problems, by = c('item_id')))
  }

 
  models = lapply(respData, fit_inter)
  
  tt = lapply(models, function(x)
  {
    A = c(1:nrow(x$ss$il))[x$ss$il$item_id %in% prop_design[prop_design[[item_property]]==props[1],]$item_id]
    B = c(1:nrow(x$ss$il))[x$ss$il$item_id %in% prop_design[prop_design[[item_property]]==props[2],]$item_id]
    SSTable(x, AB = list(A,B), model = model)
  })
  
  
  maxA = nrow(tt[[1]]$tbl)-1
  maxB = ncol(tt[[1]]$tbl)-1
  sg = data.frame(k=0:(maxA+maxB))
  
  default.args = list(asp=1, main="Profile plot", xlab=props[1], 
                      ylab=props[2],xlim=c(0,maxA),ylim=c(0,maxB))
  do.call(graphics::plot, 
          merge_arglists(user.args, 
                         default=default.args,
                         override=list(x=c(0,maxA), y=c(0,maxB),
                         xaxs="i", type="n")))
  
  # The timolines
  k = maxA + maxB
  sg$y0 = pmin(maxB,sg$k)
  sg$x0 = sg$k - sg$y0
  sg$x1 = pmin(maxA,sg$k)
  sg$y1 = sg$k - sg$x1
  graphics::segments(sg$x0, sg$y0, sg$x1, sg$y1, col="gray")
  
  graphics::text(0:maxA,0,0:maxA,cex=.6,col="lightgray")
  graphics::text(maxA,1:maxB,(maxA+1:maxB),cex=.6,col="lightgray")
  

  for (i in seq_along(tt)) {
    ta = tt[[i]]$tbl
    y = tibble(
      value=as.vector(ta),
      Var1=as.integer(gl(nrow(ta),1,nrow(ta)*ncol(ta))),
      Var2=as.integer(gl(ncol(ta),nrow(ta),nrow(ta)*ncol(ta))),
      v = as.integer(gl(nrow(ta),1,nrow(ta)*ncol(ta))) + 
        as.integer(gl(ncol(ta),nrow(ta),nrow(ta)*ncol(ta)))
    )
    stp = y %>% 
      group_by(.data$v) %>%
      do(.[which.max(.$value),]-1)
    graphics::lines(stp$Var1, stp$Var2, col=i+1, lw=2)
  }
  
  graphics::legend("topleft", 
                   legend=names(tt), 
                   lty=1, col=1+(1:length(tt)),
                   cex=.7, 
                   box.lty=0)
  graphics::box()
}

#' Exploratory DIF test
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param covariate Person covariate name for subsetting
#' @return An object of class \code{DIF_stats} holding statistics for
#' overall-DIF and a matrix of statistics for DIF in the relative position of
#' item-category parameters in the regular parameterization used e.g., by OPLM.
#' @details 
#' Tests for DIF as described in Bechger and Maris (2014; A Statistical Test 
#' for Differential Item Pair Functioning. Psychometrica). Supplements the 
#' confirmatory approach of the profile plot
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' dd = DIF(db,covariate="gender")
#' print(dd)
#' plot(dd)
#' }
#' 
DIF = function(dataSrc, covariate, predicate=NULL) 
{
  ## 1. Interpret input.. much like beginning of profile plot
  # Check whether there are 2 groups 
  # Check whether predicate (if not NULL) does not leave us without data
  if(!inherits(dataSrc,'data.frame'))
  {
    covariate = tolower(covariate)
  }
  
  columns = c('person_id','item_id','item_score', covariate)
  
  qtpredicate=eval(substitute(quote(predicate)))
  respData = get_responses_(dataSrc, qtpredicate, columns = columns, env = caller_env()) 
  if(nrow(respData) == 0) stop('no data to analyse')

  # split into list
  respData = by(respData, respData[[covariate]], identity)
  
  if(length(respData)!=2)
    stop('The covariate needs to have two unique values in your data to calculate DIF')
  
  common_scores = Reduce(dplyr::intersect, lapply(respData, function(x) x[,c('item_id','item_score')]))
  problems = Reduce(dplyr::union,
                    lapply(respData,
                           function(data) dplyr::setdiff(data[,c('item_id','item_score')],common_scores)))
  
  if(nrow(problems) > 0)
  {
    problems = tibble(item_id = unique(problems$item_id))
    warning(paste('the following items do not have the same score categories over both covariates and',
                  'have been removed from the analysis:', paste0(problems$item_id,collapse=', ')))
    respData = lapply(respData, function(x) x %>% anti_join(problems, by = 'item_id'))
  }
  
  ## 2. Estimate models with fit_enorm using CML
  models = lapply(respData, fit_enorm)
  
  ## 3. Make sure parameters pertain to same items-responses in the same order
  # This should always be correct, since fit_enorm orders on items and scores and
  # I made sure in the prelim that both models have the same items and score catregories
  
  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = OverallDIF_ (models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
                          models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)

  D = PairDIF_(models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
                models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)
  
  
  ## 5. Report D and DIF_stats and inputs
  ou = list(DIF_overall = DIF_stats, DIF_pair = D, 
            groups = names(respData), items = unique(common_scores$item_id))
  class(ou) = append('DIF_stats', class(ou))
  return(ou)
}


print.DIF_stats <- function(x, ...)
{
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  tmp = specify_decimal(x$DIF_overall$p,3)
  if (tmp=="0.000") tmp="< 0.0006"
  print(paste0("Test for DIF:"," Chi-square = ", as.character(round(x$DIF_overall$stat, digits=3)),
               ", df = ", 
               as.character(x$DIF_overall$df),
               ", p = ", tmp))  
}

plot.DIF_stats = function(x, ...)
{
  D = x #had to make the parameter name x instead of D for CRAN rules
  x=D$DIF_pair
    yLabels = rownames(x)
    xLabels = colnames(x)
    min_=min(x)
    max_=max(x)
    default.args = list(main = paste(D$groups, collapse = ' vs. '),
                        axes=FALSE, zlim=c(min_,max_),xlab='',ylab='')
    
    graphics::layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
    

    tmp = rainbow(256)[1:128]
    ColorRamp=c(tmp, tmp[128:1])
    ColorLevels <- seq(min(x), max(x), length=length(ColorRamp))
    
    # Reverse Y axis
    reverse <- nrow(x) : 1
    yLabels <- yLabels[reverse]
    x <- x[reverse,]
    
    # Data Map
    oldpar = par(mar = c(6,8,2.5,2))
    do.call(image,
            merge_arglists(list(...),
                           override = list(x = 1:length(xLabels), y = 1:length(yLabels), z=t(x),
                                           col=ColorRamp),
                           default = default.args))


    axis(1, at=1:length(xLabels), labels=xLabels, las= 3, cex.axis=0.6)
    axis(2, at=1:length(yLabels), labels=yLabels, las=1,
         cex.axis=0.6)
    
    #Color Scale
    par(mar=c(6,3,2,2))
    image(1, ColorLevels,
            matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
            col=ColorRamp,
            xlab="",ylab="",
            xaxt="n")

    graphics::layout(1)
    par(oldpar)
  
}


############################
#' Simple test-item analysis
#'
#' Show simple Classical Test Analysis statistics
#' at item and test level
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param type How to present the item level statistics: \code{raw} for each test booklet 
#' separately, \code{averaged} averaged over the test booklet in which the item is included,
#' with the number of persons as weights, or {compared}, in which case the pvalues, 
#' correlations with the sum score (rit), and correlations with the rest score (rit) are 
#' shown in separate tables and compared across booklets
#' @return A list containing \item{testStats}{a data frame of statistics at test level}, 
#' and \item{itemStats}{a data frame of statistics at item level}.
#'
tia_tables <- function(dataSrc, predicate = NULL, type=c('raw','averaged','compared')) {
  type = match.arg(type)

  qtpredicate = eval(substitute(quote(predicate)))
  x = get_sumscores(dataSrc, qtpredicate, env=caller_env())
  if(nrow(x) == 0) stop('no data to analyse')

  maxScores = x %>%
    group_by(.data$item_id) %>%
    summarise(maxScore=max(.data$item_score))

  ti = x %>% ungroup() %>% group_by(.data$booklet_id,.data$item_id) %>%  
    summarise(meanScore=mean(.data$item_score),
              sdScore=sd(.data$item_score),
              rit=cor(.data$item_score, .data$sumScore),
              rir=cor(.data$item_score, .data$sumScore - .data$item_score),
              n=n()) %>% 
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
                                    Pval=weighted.mean(.data$pvalue, w=.data$n, na.rm=TRUE),
                                    Rit=weighted.mean(.data$rit, w=.data$n, na.rm=TRUE),
                                    Rir=weighted.mean(.data$rir, w=.data$n, na.rm=TRUE),
                                    N=sum(.data$n, na.rm=TRUE),
                                    mnScore=weighted.mean(.data$meanScore, w=.data$n, na.rm=TRUE),
                                    sdScore=weighted.mean(.data$sdScore, w=.data$n, na.rm=TRUE))
                     },
                     compared={
                       list(
                         pvalue = ti %>% select(.data$booklet_id,.data$item_id,.data$pvalue) %>% spread_(key_col='booklet',value_col='pvalue'),
                         rit =    ti %>% select(.data$booklet_id,.data$item_id,.data$rit) %>% spread_(key_col='booklet',value_col='rit'),
                         rir =    ti %>% select(.data$booklet_id,.data$item_id,.data$rir) %>% spread_(key_col='booklet',value_col='rir')
                       )
                     }
  )
  testStats = ti %>% ungroup() %>%
    group_by(.data$booklet_id) %>% 
    filter(complete.cases(.data$sdScore)) %>%
    summarise(nItems=n(),
              alpha=.data$nItems/(.data$nItems-1)*(1-sum(.data$sdScore^2) / sum(.data$rit * .data$sdScore)^2 ),
              meanP=mean(.data$pvalue),
              meanRit=mean(.data$rit),
              meanRir=mean(.data$rir),
              maxTestScore=sum(.data$maxScore),
              N=max(.data$n))
     
  list(itemStats=itemStats, testStats=testStats)
}

  
###########################################################  
#' Interactive test-item analysis
#'
#' Opens up a shiny application to do interactive item-test
#' analysis on the database
#'
#'
#' @param db A handle to the database, i.e. the output of \code{create_new_project}
#' or \code{open_project}
#'
iTIA <- function(db) {
  
  theTIA = tia_tables(db, type='averaged')
  ourItems = theTIA$itemStats$item_id
  ourBooklets = theTIA$testStats$booklet_id
    
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Booklets", tabName = "booklets"),
      menuItem("Items", tabName = "items")
    )
  )
  
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "booklets", dataTableOutput("booklets"),
              lapply(ourBooklets,
                     function(i){
                       bsModal(paste0("myModal", i),
                               # "Local independence (left), Interaction model (right)",
                               "Models", 
                               paste0("btn", i),
                               size = "large",
                               plotOutput(paste0("plot", i)))
                     })
      ),
      tabItem(tabName = "items", dataTableOutput("items"),
              lapply(ourItems, function(i){
                bsModal(paste0("myModal", i),
                        "Distractor plot",
                        paste0("btn", i),
                        size = "large",
                        plotOutput(paste0("plot", i)))
              })
      )
    )
  )
  
  # Put them together into a dashboardPage
  ui = dashboardPage(
    dashboardHeader(title = "Interactive Test-Item Analysis"),
    sidebar,
    body
  )
  
  server <- function(input, output, session) {
    tia = reactive({
      #result = tia_tables(db, type='averaged')
      #result
      theTIA
    })
    
    lapply(ourItems, function(i){
      output[[paste0("plot", i)]] =
        renderPlot(distractor_plot(db,i,nr=3,nc=3))
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
    
    lapply(ourBooklets, function(i){
      mo = fit_inter(db, booklet_id==i)
      output[[paste0("plot", i)]] = renderPlot(
        plot(mo, overlay=TRUE, nc=2)
      )
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
    
    output$booklets = DT::renderDataTable({
      B=tia()[[2]]
      B$alpha=round(B$alpha,3)
      B$meanP=round(B$meanP,3)
      B$meanRit=round(B$meanRit,3)
      B$meanRir=round(B$meanRir,3)
      Plots = shinierInput(actionButton, "btn", ourBooklets, label = "Show")
      B = cbind(B, Plots)
      B
    }, extensions = 'Buttons', options=list(dom = 'Bfrtip',
                                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
                                            orderClasses = TRUE,
                                            preDrawCallback = DT::JS("function() {
                                                                     Shiny.unbindAll(this.api().table().node()); }"),
                                            drawCallback = DT::JS("function() {
                                                                  Shiny.bindAll(this.api().table().node()); } "),
                                            lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
                                            autoWidth = TRUE, scrollX=TRUE), escape=FALSE)
    
    
    output$items = DT::renderDataTable({
      A=tia()[[1]]
      A$Pval=round(A$Pval,3)
      A$Rit=round(A$Rit,3)
      A$Rir=round(A$Rir,3)
      A$mnScore=round(A$mnScore,2)
      A$sdScore=round(A$sdScore,2)
      Plots <- shinierInput(actionButton, "btn", ourItems, label = "Show")
      A = cbind(A,Plots)
      A
    }, extensions = 'Buttons', options = list(dom = 'Bfrtip',
                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
                                              orderClasses = TRUE,
                                              pageLength = 25,
                                              preDrawCallback = DT::JS("function() {
                                                                       Shiny.unbindAll(this.api().table().node()); }"),
                                              drawCallback = DT::JS("function() {
                                                                    Shiny.bindAll(this.api().table().node()); } "),
                                              lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
                                              autoWidth = TRUE,
                                              scrollX=TRUE
                                              ), escape=FALSE)
    } # end of server
  shinyApp(ui, server)
}


#############################
#' Interactive model display
#'
#' Opens up a shiny application with item statistics and interactive
#' plots for the Rasch and Interaction models
#'
#'
#' @param db A handle to the database, i.e. the output of \code{create_new_project}
#' or \code{open_project}
#' @param booklet ID of the booklet that will be shown
#' 
iModels <- function(db, booklet) {
  models = fit_inter(db, booklet_id==booklet)
  tia = tia_tables(db, booklet_id==booklet, type='raw')$itemStats 
  sidebar = dashboardSidebar(
    checkboxInput("show", "Show data", value = FALSE, width = NULL),
    checkboxInput("summate", "As scores", value = TRUE, width = NULL)
  )
  body = dashboardBody(
    dataTableOutput("items"),
    lapply(tia$item_id,
           function(i){
             bsModal(paste0("myModal", i),
                     "Models",
                     paste0("btn", i),
                     size = "large",
                     plotOutput(paste0("plot", i)))
           })
  )
  # Put them together into a dashboardPage
  ui = dashboardPage(
    dashboardHeader(title = "Interactive Models"),
    sidebar,
    body
  )
  
  server <- function(input, output, session) {
    output$items = DT::renderDataTable({
      A=tia[tia$booklet_id==booklet,c("item_id","meanScore","maxScore","pvalue","rit","rir","n")]
      A$pvalue=round(A$pvalue,3)
      A$rit=round(A$rit,3)
      A$rir=round(A$rir,3)
      A$meanScore=round(A$meanScore,2)
      Plots = shinierInput(actionButton,  "btn", tia$item_id, label = "Show")
      A = cbind(A,Plots)
      A
    },extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
      orderClasses = TRUE,
      pageLength = 25,
      preDrawCallback = DT::JS("function() {
                               Shiny.unbindAll(this.api().table().node()); }"),
      drawCallback = DT::JS("function() {
                            Shiny.bindAll(this.api().table().node()); } "),
      lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
      autoWidth = TRUE,
      scrollX=TRUE
      ), escape=FALSE)
    lapply(tia$item_id, function(i){
      output[[paste0("plot", i)]] =
        renderPlot(plot(models, item=i,
                        summate=input$summate,
                        overlay=FALSE,
                        nc=1, nr=1, curtains=10,
                        show.observed=input$show
        ))
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
      }
  shinyApp(ui, server)
}





##########################################
#' Fit the extended nominal response model
#'
#' Fits the extended nominal response model by Bayesian sampling or CML
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nIterations Number of Gibbs samples when method is Bayes, max. number of iterations 
#' when method is CML
#' @return Depends on the estimation method
#'
fit_enorm <- function(dataSrc, predicate = NULL, method=c("CML", "Bayes"), nIterations=500) {
  
  method <- match.arg(method)
  qtpredicate = eval(substitute(quote(predicate)))

  x = get_sumscores(dataSrc, qtpredicate, env=caller_env())
  if(nrow(x) == 0) stop('no data to analyse')
  
  design = get_design_from_sumscores(x) %>% select(-.data$n_persons)
  if(length(unique(design$booklet_id)) > 1)
  {
    im = as.matrix(table(design$item_id, design$booklet_id))
    wm = crossprod(im, im)
    diag(wm) = 0
    if(!design_is_connected(list(im=im, wm=wm))) stop('Your design is not connected')  
  }
  
  # new dplyr requires ungroup
  itm_max = x %>% group_by(.data$item_id) %>% summarise(maxScore=max(.data$item_score)) 
  #design =  x %>% group_by(.data$booklet_id, .data$item_id) %>% slice(1)  %>% select(.data$booklet_id,.data$item_id)
  
  ssBIS = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$item_score) %>% 
    summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>% 
    ungroup()
  
  plt = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$sumScore) %>% 
    summarise(meanScore=mean(.data$item_score), N=n()) %>% 
    ungroup()

  maxScores = itm_max %>%
    inner_join(design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    summarise(maxTotScore = sum(.data$maxScore))
  
  allScores = maxScores %>% group_by(.data$booklet_id) %>%
    do(tibble(sumScore=0:.$maxTotScore))
  
  stb = plt %>%
    select(.data$booklet_id, .data$sumScore, .data$N) %>%
    distinct() %>%
    right_join(allScores, by=c('booklet_id','sumScore')) %>%
    do({
      .$N[is.na(.$N)]=0
      .$booklet_id[is.na(.$booklet_id)] = .$booklet_id[!is.na(.$booklet_id)][1]
      as.data.frame(.)
    }) %>%
    select(.data$booklet_id, .data$sumScore, .data$N) %>%
    arrange(.data$booklet_id, .data$sumScore)
  
  ssIS = ssBIS %>% 
    group_by(.data$item_id, .data$item_score) %>%
    summarise(sufI = sum(.data$sufI)) %>%
    arrange(.data$item_id, .data$item_score)
  
  ssI  = ssIS %>% 
    summarise(nCat = n()) %>% 
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1,last = cumsum(.data$nCat)) %>%
    arrange(.data$item_id)

  m = x  %>% 
    group_by(.data$booklet_id)  %>% 
    summarise(m = n_distinct(.data$person_id))
  
  a = ssIS$item_score
  b = exp(runif(nrow(ssIS), -1, 1))

  bkl = lapply(m$booklet_id, function(x) {
    itInBk = design$item_id[design$booklet_id==x] 
    items = ssI$item_id[ssI$item_id %in% itInBk]
    m = m$m[m$booklet_id==x]
    first = ssI$first[match(items, ssI$item_id)]
    last =  ssI$last[match(items, ssI$item_id)]
    list(booklet=x, items=items, m=m, first=first, last=last, 
         score=stb$sumScore[stb$booklet_id==x], 
         scoretab=stb$N[stb$booklet_id==x], 
         lambda=rep(1,sum(a[last])+1))
  })
  
  names(bkl) = m$booklet_id
  
  itemList = lapply(ssI$item_id, function(x) design$booklet_id[design$item_id==x])
  
  itemListInt = lapply(itemList,function(x)match(x,m$booklet_id))

  it_sc_lab=paste0(ssIS$item_id[-ssI$first], "_",ssIS$item_score[-ssI$first])
  if (method=="CML"){
    result = try(calibrateCML(booklet=bkl, sufI=ssIS$sufI, a=ssIS$item_score, first=ssI$first, last=ssI$last, nIter=nIterations))
    names(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score)
    row.names(result$beta.cml)=it_sc_lab
    rownames(result$acov.cml)=it_sc_lab
    colnames(result$acov.cml)=it_sc_lab
  } else {
    result = try(calibrate(itemList=itemListInt, booklet=bkl, sufI=ssIS$sufI, b=b, a=a, first=ssI$first, last=ssI$last, nIter=nIterations))
    colnames(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score) 
    colnames(result$beta.cml)=it_sc_lab
  }
  if (inherits(result, "try-error")) result=NULL

  outpt = list(est=result, inputs=list(bkList=bkl, ssIS=ssIS %>% ungroup(), ssI=ssI, stb=stb, method=method), xpr=as.character(qtpredicate))
  class(outpt) = append('prms', class(outpt))
  outpt
}



##########################################
#' Draw plausible values
#'
#' Draws plausible values based on sum scores and a fitted
#' ENORM model
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param nPV Number of plausible values to draw per person.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler (this is 
#' recognised automatically), the number of the random draw (iteration) to use 
#' in generating the PV. If NULL, all draws will be averaged. If outside range,
#' the last iteration will be used.   
#' @return Depends on the estimation method
#' 
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'    covariates=list(gender="<unknown>"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' f=fit_enorm(db)
#' par(mfrow=c(1,2))
#' pv_M=plausible_values(db,f,(mode=="Do")&(gender=="Male"))
#' pv_F=plausible_values(db,f,(mode=="Do")&(gender=="Female"))
#' plot(ecdf(pv_M$PV), 
#'    main="Do: males versus females", xlab="Ability", col="red")
#' lines(ecdf(pv_F$PV), col="green")
#' legend(-2.2,0.9, c("female", "male") , 
#'    lty=1, col=c('green', 'red'), bty='n', cex=.75)
#'
#' pv_M=plausible_values(db,f,(mode=="Want")&(gender=="Male"))
#' pv_F=plausible_values(db,f,(mode=="Want")&(gender=="Female"))
#'
#' plot(ecdf(pv_M$PV), 
#'    main="Want: males versus females", xlab=" Ability", col="red")
#' lines(ecdf(pv_F$PV),col="green")
#' legend(-2.2,0.9, c("female", "male") , 
#'    lty=1, col=c('green', 'red'), bty='n', cex=.75)
#' }
#' 
plausible_values <- function(dataSrc, parms, predicate=NULL, nPV=1, use_draw=NULL)
{
  qtpredicate=eval(substitute(quote(predicate)))
  plausible_values_(dataSrc, parms, qtpredicate=qtpredicate, nPV=nPV, use_draw=use_draw, env=caller_env())
}
  
plausible_values_ <- function(dataSrc, parms, qtpredicate=NULL, nPV=1, use_draw=NULL, env=NULL)
{
  if(is.null(env)) env = caller_env()
  x = get_sumscores(dataSrc, qtpredicate, env=caller_env())
  if(nrow(x) == 0) stop('no data to analyse')
  # we make the design and join with the params
  design =  x %>% 
    group_by(.data$booklet_id, .data$item_id) %>% 
    slice(1) %>% 
    select(.data$booklet_id, .data$item_id) %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  x = x %>% group_by(.data$booklet_id, .data$person_id) %>% slice(1)
  
  if(parms$input$method=='CML') {
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
    #bkl = parms$inputs$bkList
  } else {
    a = parms$est$a
    #bkl = parms$inputs$bkList
    if(is.null(use_draw)) {
      b = colMeans(parms$est$b)  
    } else {
      if (use_draw %in% 1:nrow(parms$est$b)) {
        b = parms$est$b[use_draw,]   
      } else {
        b = parms$est$b[nrow(parms$est$b),]
      }
    }   
  } # extracted parms
  y = x %>% ungroup() %>%
    group_by(.data$booklet_id) %>%
    do({
      bkID = .$booklet_id[1]
      first=design[design$booklet_id==bkID,]$first
      last=design[design$booklet_id==bkID,]$last
      PV = pv(b,a,first,last, .$sumScore,nPV) 
      data.frame(.$person_id, .$sumScore, as.data.frame(PV))
    })
  names(y) = c('booklet_id','person_id','sumScore',paste0('PV',1:nPV))
  y
}




##########################################
#' Estimate abilities
#'
#' Computes MLE of ability and optionally attaches them to person data by sum score
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method If ability will be estimated with maximum likelihood (MLE). Otherwise, 
#' we produce a Bayesian expected a posteriori (EAP) estimate.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler (this is 
#' recognised automatically), the number of the random draw (iteration) to use 
#' in generating the PV. If NULL, all draws will be averaged. If outside range,
#' the last iteration will be used.
#' @param person_level If TRUE, return results per person, otherwise just
#' the score transformation table.   
#' @return Depends on \code{person_level}
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' f = fit_enorm(db)
#' aa = ability(db,f,method="MLE",person_level = FALSE)
#' bb = ability(db,f,method="EAP",person_level = FALSE)
#' plot(bb$sumScore, bb$theta, xlab="test-score", ylab="ability est.", pch=19, cex=0.7)
#' points(aa$sumScore, aa$theta, col="red", pch=19, cex=0.7)
#' legend("topleft", legend = c("EAP", "MLE"), bty = "n",
#'     lwd = 1, cex = 0.7, col = c("black", "red"), lty=c(0,0), pch = c(19,19))
#' }
#' 
ability <- function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP"), use_draw=NULL, person_level=TRUE){

  method <- match.arg(method)
  qtpredicate=eval(substitute(quote(predicate)))
  x = get_sumscores(dataSrc, qtpredicate, env=caller_env())
  if(nrow(x) == 0) stop('no data to analyse')
  
  # we make the design and join with the params
  design =  x %>% 
    group_by(.data$booklet_id, .data$item_id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(.data$booklet_id, .data$item_id) %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id,.data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  x = x %>% group_by(.data$booklet_id, .data$person_id) %>% slice(1)
  
  if(parms$input$method=="CML"){
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
    #bkl = parms$inputs$bkList
  } else {
    a = parms$est$a
    #bkl = parms$inputs$bkList
    if(is.null(use_draw)) {
      b = colMeans(parms$est$b)  
    } else {
      if (use_draw %in% 1:nrow(parms$est$b)) {
        b = parms$est$b[use_draw,]   
      } else {
        b = parms$est$b[nrow(parms$est$b),]
      }
    }   
  } # extracted parms
  # is the maxscore always sum(nCat-1), Ivaylo?
  #booklets = design %>% group_by(.data$booklet_id) %>% summarise(maxScore=sum(.data$nCat-1))
  #ou = design %>% group_by(.data$booklet_id) %>% summarise(maxScore=sum(.data$nCat-1))
  
  if (method=="MLE")
  {
    mx_scores = parms$inputs$ssIS %>%  group_by(.data$item_id) %>% summarise(item_max=max(.data$item_score))
    ou = design %>% inner_join(mx_scores, by='item_id') %>% group_by(.data$booklet_id) %>% summarise(maxScore=sum(.data$item_max))
    ou = merge(ou,data.frame(sumScore=0:max(ou$maxScore)))
    ou = ou[ou$sumScore <= ou$maxScore,c('booklet_id','sumScore')] %>% arrange(.data$booklet_id, .data$sumScore)
    ou$theta = unlist(map(sort(unique(ou$booklet_id)), ~theta_MLE(b,a,design[design$booklet_id==.x,]$first,
                                                             design[design$booklet_id==.x,]$last)))
  }else
  {
    mx_scores = parms$inputs$ssIS %>%  group_by(.data$item_id) %>% summarise(item_max=max(.data$item_score))
    ou = design %>% inner_join(mx_scores, by='item_id') %>% group_by(.data$booklet_id) %>% summarise(maxScore=sum(.data$item_max))
    ou = merge(ou,data.frame(sumScore=0:max(ou$maxScore)))
    ou = ou[ou$sumScore <= ou$maxScore,c('booklet_id','sumScore')] %>% arrange(.data$booklet_id, .data$sumScore)
    ou$theta = unlist(map(sort(unique(ou$booklet_id)), ~theta_EAP(b,a,design[design$booklet_id==.x,]$first,
                                                                   design[design$booklet_id==.x,]$last,ou$sumScore[ou$booklet_id==.x])))

  }
  if (!person_level) return(ou) else{
    return(x %>% inner_join(ou, by = c("booklet_id", "sumScore")) %>% 
           select(.data$booklet_id, .data$person_id, .data$sumScore, .data$theta))
  }
}


