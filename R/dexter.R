
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

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
#' @param covariates An optional list of person covariates. Names should correspond to covariates intended to be used in the project.
#' Values are used as default (missing) values for these covariates. The datatype will also be inferred from the values.
#' Known covariates will be imported (if supplied) in \code{\link{add_booklet}}. 
#' @return If the scoring rules pass a sanity check, a handle to the data base.
#' Otherwise, a data frame listing the problems found, with 4 columns:
#' item_id: id of the problematic item
#' less_than_two_scores: if TRUE, the item has only one distinct score
#' duplicated_responses: if TRUE, the item contains two or more identical response categories
#' min_score_not_zero: if TRUE, the minimum score of the item was not 0
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
#'\donttest{
#' head(verbAggrRules)
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#'   covariates=list(gender="<unknown>"))
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
  
  if(any(is.na(rules$item_id)) || any(is.na(rules$item_score)))
    stop("The item_id and item_score columns may not contain NA values")
  
  sanity = rules %>%
    group_by(.data$item_id) %>%
    summarise(less_than_two_scores = length(unique(.data$item_score))<2,
              duplicated_responses = any(duplicated(.data$response)),
              min_score_not_zero = min(.data$item_score)>0) %>%
    filter(.data$less_than_two_scores | .data$duplicated_responses | .data$min_score_not_zero)

  if (nrow(sanity)>0) {
    message("There were problems with your scoring rules.\nCheck the output below for possible reasons.\n")
    print(as.data.frame(sanity))
    stop('Scoring rules are not valid')
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
    },on_error = function(e){dbDisconnect(db);stop(e)})
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
  
  if (!file.exists(db_name)) stop("There is no such file")
    
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
                     'to attempt to convert your database to the newest version of Dexter.'))
      } else 
      { 
        message('Converting your database to the new Dexter version.')
        convert_old_db(db)
      }
    } else
    {
      dbDisconnect(db)
      stop('Sorry, this does not appear to be a Dexter database.')
    }
  }
  return(db)
}

#' Close a project
#'
#' This is just an alias for \code{DBI::dbDisconnect(db)}, included for completeness
#'
#' @param db a Dexter project db handle
#'
close_project = function(db) dbDisconnect(db) 


##################################
#' Derive scoring rules from keys
#'
#' For multiple choice items that will be scored as 0/1, derive the
#' scoring rules from the keys to the correct responses
#'
#'
#' @param keys  A data frame containing columns \code{item_id}, \code{nOptions}, and
#' \code{key} (the spelling is important). See details.
#' @param include_NA_rule whether to add an option 'NA' (which is scored 0) to each item
#' @return A data frame that can be used as input to \code{start_new_project}
#' @details
#' This function might be useful in setting up the scoring rules when all items
#' are multiple-choice and scored as 0/1. (Hint: Because the order in which the
#' scoring rules is not important, one can use the function to generate rules for
#' many MC items and then append per hand the rules for a few complex items.)
#'
#' The input data frame must contain the exact name of each item, the number
#' of options, and the key. If the keys are all integers, it will be assumed that
#' responses are coded as 1 through nOptions. If they are all uppercase letters,
#' it is assumed that responses are coded as A,B,C,... All other cases result
#' in an error.
#'
keys_to_rules <- function(keys, include_NA_rule = FALSE) {
  # for backward compatibility we rename
  colnames(keys)[colnames(keys)=='item'] = 'item_id'
  
  keys = mutate_if(keys, is.factor,as.character)
  
  if (is.numeric(keys$key)) ABC=FALSE else {
    if (all(keys$key %in% LETTERS)) ABC=TRUE
    else stop("You have inadmissible keys")
  }
  if (ABC) {
    m = match(keys$key, LETTERS)
    if (any(m>keys$nOptions)) stop("You have out-of-range keys")
    r = keys %>% 
      group_by(.data$item_id) %>% 
      do({
      y = tibble(response=LETTERS[1:.$nOptions[1]], item_score=0)
      y$item_score[match(.$key[1],LETTERS)] = 1
      y
    })
  } else {
    if (any(keys$key>keys$nOptions)) stop("You have out-of-range keys")
    r = keys %>% group_by(.data$item_id) %>% do({
      y = tibble(response=1:.$nOptions[1], item_score=0)
      y$item_score[.$key[1]] = 1
      y
    })
  }
  r = ungroup(r)
  
  if(include_NA_rule){
    r = bind_rows(r, tibble(item_id=unique(r$item_id), response='NA', item_score=0))
  }
  r$item_score = as.integer(r$item_score)
  r
}



#' Add or modify scoring rules
#' 
#' Having to alter or add a scoring rule is occasionally necessary, e.g. in case of a key error. 
#' This function offers the possibility to do so and also allows you to add new items to your project
#' 
#' @param db handle to a Dexter project database
#' @param rules A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
#' The order is not important but spelling is. Any other columns will be ignored. See details
#' @return If the scoring rules pass a sanity check, a small summary of changes is printed and nothing is returned
#' Otherwise this function returns a data frame listing the problems found, with 4 columns:
#' item_id: id of the problematic item
#' less_than_two_scores: if TRUE, the item has only one distinct score
#' duplicated_responses: if TRUE, the item contains two or more identical response categories
#' min_score_not_zero: if TRUE, the minimum score of the item was not 0
#' Please note that the sanity check is done for the items you change

#' @details 
#' The rules should contain all rules that you want to change or add. This means that in case of a key error
#' in a single multiple choice question, you typically have to change two rules.
#' @examples 
#'\donttest{
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
  
  if(any(is.na(rules$item_id)) || any(is.na(rules$item_score)))
    stop("The item_id and item_score columns may not contain NA values")
  
  rules = rules[, c("item_id", "response", "item_score")]
  rules$response = as.character(rules$response)
  rules$response[is.na(rules$response)] = 'NA'
  rules$item_id = as.character(rules$item_id)
  
  # check if same options are supplied multiple times, this is the only check not done in conjuntion
  # with the existing rules so it has to be gotten out of the way
  items_with_duplicate_options = unique(rules[duplicated(select(rules,.data$item_id,.data$response)),]$item_id)

  existing_rules = dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxScoring_rules;')
  # remove the no-ops
  rules = dplyr::setdiff(rules, existing_rules)
  
  existing_opts = existing_rules %>% select(-.data$item_score)
  
  # new items or responses
  new_rules = rules %>% anti_join(existing_opts, by=c('item_id','response'))
  
  # existing items or responses but with new scores
  amended_rules = rules %>% inner_join(existing_opts, by=c('item_id','response'))
  
  # to judge the validity of the new rules, we have to look at them in combination
  # with the rules in the db that will not be changed
  sanity = rules %>% 
    dplyr::union(existing_rules %>% 
                   anti_join(rules, by=c('item_id','response'))
                 ) %>%
    group_by(.data$item_id) %>%
    summarise(less_than_two_scores = length(unique(.data$item_score))<2,
              duplicated_responses = any(duplicated(.data$response)),
              min_score_not_zero = min(.data$item_score)>0) 
  
  if(any(sanity$item_id %in% items_with_duplicate_options))
    sanity[sanity$item_id %in% items_with_duplicate_options,]$duplicated_responses = TRUE
  
  sanity = filter(sanity, .data$less_than_two_scores | .data$duplicated_responses | .data$min_score_not_zero)
  
  
  if (nrow(sanity)) {
    message("There were problems with your scoring rules.\nCheck the output below for possible reasons.\n")
    print(as.data.frame(sanity))
    stop('Scoring rules are not valid')
  }    
  
  dbTransaction(db,
  {
    if(nrow(new_rules)>0)
    {
      new_items = setdiff(new_rules$item_id, dbGetQuery(db, 'SELECT item_id FROM dxItems;')$item_id)
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
  message(paste0('rules_changed: ', nrow(amended_rules), '\nrules_added: ', nrow(new_rules)))
}


#' Add a booklet to a project
#'
#' Adds item response data for a test form (a.k.a. booklet)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param x A data frame containing the responses and, optionally,
#' person covariates. The data.frame should have one row per respondent and the column names should 
#' correspond to the item_id's in the rules or the names of the covariates. See details.
#' @param booklet_id A (short) string identifying the test form (booklet)
#' @param auto_add_unknown_rules  If FALSE, an error will be generated if 
#' some of the responses do not appear in the scoring rules. Default is TRUE.
#' @return A list of: \item{items}{The names of the columns in \code{x} that were
#' treated as items}
#' \item{covariates}{The names of the columns in \code{x} that were
#' treated as person covariates}
#' \item{not_listed}{A data frame of all responses that will be treated as missing}
#' @details It is a common practice to keep respons data in tables where each row 
#' contains the responses from a single person. This function is provided to input
#' data in that form, one booklet at a time. The starting point is a data frame.
#' Getting the data frame into R is left to the user. We have found packages
#' \code{readxl} to be very good at reading Excel sheets, and \code{haven} quite
#' efficient with SPSS files.
#'
#' If the dataframe \code{x} contains a variable named \code{person_id} this variable 
#' will be used to identify unique persons. It is assumed that a single person will only 
#' make a single booklet once, otherwise an error will be generated. 
#' 
#' If a person_id is not supplied, dexter will generate unique person_id's for each row of data.  
#'
#' Any column whose name has an exact match in the scoring rules inputted with
#' function \code{start_new_project} will be treated as an item; any column whose name has an 
#' exact match in the covariates will be treated as covariate. If a name matches both
#' a covariate and an item, the item takes precedence. Variables other than items, covariates 
#' and person_id will be ignored.
#' 
#' If \code{auto_add_unknown_rules=TRUE}, any responses to an item that do not have an 
#' exact match in the scoring rules will be automatically given the lowest score of 0. 
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
#' 
#' close_project(db)
#' }
#' 
add_booklet <- function(db, x, booklet_id, auto_add_unknown_rules = TRUE) {
  x = x %>% mutate_if(is.factor, as.character) 
  
  covariates = intersect(dbListFields(db, 'dxPersons'), tolower(names(x)))
  covariates = covariates[covariates != 'person_id']
  
  design = tibble(booklet_id = booklet_id, item_id = names(x), col_order=c(1:ncol(x))) %>%
    inner_join(dbGetQuery(db, "SELECT item_id FROM dxItems;"), by='item_id') %>%
    mutate(item_position = dense_rank(.data$col_order)) %>%
    select(-.data$col_order)

  if(nrow(design) == 0) stop('None of the column names in x correspond to known items in your project')

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
      auto_add_unknown_rules = auto_add_unknown_rules,
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
#' @param overwrite Whether existing item properties should be overwritten (default=TRUE)
#' @param default_values a list where the names are item_properties and the values are defaults.
#' The defaults will be used as the value for an item for which the property is not specified, 
#' for example when you add new items using \code{\link{touch_rules}}.
#' Default_values for an item property will only be processed
#' the first time you define an item property.
#' @return nothing
#' @details When entering response data in the form of a rectangular person x item
#' table, it is easy to provide person properties but practically impossible
#' to provide item properties. This function provides a possibility to do so.
#' The order of the rows and columns in the data frame is not important but
#' there must be a column called \code{item_id} containing the item id's
#' exactly as entered before.
#' 
#' Note that is is not possible to add new items with this function, 
#' use \code{\link{touch_rules}} if you want to add new items to your project.
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
#' get_item_properties(db) 
#' 
#' close_project(db)
#' }
#'
add_item_properties <- function(db, item_properties, overwrite=TRUE, default_values=list()) {
  item_properties = item_properties %>%
    mutate_if(is.factor, as.character) 
  
  colnames(item_properties)[colnames(item_properties)=='item'] = 'item_id'
  colnames(item_properties) = dbValid_colnames(colnames(item_properties))
  
  names(default_values) = dbValid_colnames(names(default_values))
  
  if(!'item_id' %in% colnames(item_properties))
  {
    stop("there was no column provided with name 'item_id'")
  }
  existing_item_properties = dbListFields(db, 'dxItems') # for convenience we include the item_id as a property
  if(!overwrite & !setequal(intersect(names(item_properties), existing_item_properties), 'item_id'))
  {
    stop('Some of the listed item properties already exist, specify overwrite=TRUE to overwrite')
  }

  #if(!setequal(dbGetQuery(db,'SELECT item_id FROM dxItems;')$item_id, item_properties$item_id))
  #{
  #  stop('properties not specified for all items.')
  #}
  dbTransaction(db, 
  {
    for(prop_name in setdiff(names(item_properties), existing_item_properties))
    {
      if(prop_name %in% names(default_values) && (is.na(default_values[[prop_name]]) && is.logical(NA)))
      {
        dbExecute(db, paste0("ALTER TABLE dxItems ADD COLUMN ",prop_name, sql_col_def(default_values[[prop_name]], TRUE, db),';'))

      } else
      {
        dbExecute(db, paste0("ALTER TABLE dxItems ADD COLUMN ",prop_name, sql_data_type(item_properties[,prop_name]),";"))
      }
    }
    pnames = names(item_properties)[names(item_properties)!='item_id']
    
    dbExecute(db,paste0('UPDATE dxItems SET ',paste0(pnames,'=:',pnames,collapse=', '),' WHERE item_id=:item_id;'),
               item_properties)
  })
  invisible(NULL)
}



#' Deprecated function names
#' 
#' All 'show_<something>' functions have been renamed to 'get_<something>'.
#' 
#' @param db handle to a Dexter project database
#' @return data.frame
#' 
#' @seealso \code{\link{get_booklets}}, \code{\link{get_design}}, \code{\link{get_item_properties}}, 
#' \code{\link{get_items}}, \code{\link{get_person_properties}}, \code{\link{get_persons}}, \code{\link{get_rules}}
show_rules = function(db)
{
  message('show_rules is deprecated, use get_rules() instead')
  get_rules(db)
}

#' @rdname show_rules
show_booklets <- function(db) {
  message('show_booklets is deprecated, use get_booklets() instead')
  get_booklets(db)
}

#' @rdname show_rules
show_item_properties <- function(db) {
  message('show_item_properties is deprecated, use get_item_properties() instead')
  get_item_properties(db)
}

#' @rdname show_rules
show_person_properties <- function(db) {
  message('show_person_properties is deprecated, use get_person_properties() instead')
  get_person_properties(db)
}

#' @rdname show_rules
show_items <- function(db){
  message('show_items is deprecated, use get_items() instead')
  get_items(db)
}

#' @rdname show_rules
show_design <- function(db){
  message('show_design is deprecated, use get_design() instead')
  get_design(db)
}

#' @rdname show_rules
show_persons <- function(db){
  message('show_persons is deprecated, use get_persons() instead')
  get_persons(db)
}


#' Get scoring rules
#' 
#' Retrieve the scoring rules currently present in the dexter project db
#' 
#' @param db handle to a Dexter project database
#' @return data.frame of scoring rules containing columns: item_id, response, item_score
#' 
get_rules = function(db)
{
  dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxScoring_rules ORDER BY item_id, response;')
}


#' Booklets entered in a project
#'
#' Retrieve information about the booklets entered in the db so far
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame showing with columns: booklet_id, n_persons and n_items.
#'
get_booklets <- function(db) {
  dbGetQuery(db,'SELECT booklet_id, n_items, n_persons FROM dxBooklet_stats ORDER BY booklet_id; ')
}


#' Item properties in a project
#'
#' Retrieve information about the item properties defined in the project (if any).
#' This will return a data.frame with one row per unique value of each item property. 
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data.frame with one row per unique value of each item property and a count of the nbr of items it applies to.
#' The columns are: item_property, value, N
#'
get_item_properties <- function(db) {
  data = dbGetQuery(db, 'SELECT * FROM dxItems;') 
  if( ncol(data) == 1)
    return(data.frame())
  
  data %>% 
    gather(key='item_property', value='value', -.data$item_id) %>%
    group_by(.data$item_property, .data$value) %>%
    summarise(N=n()) %>%
	  ungroup()
}



#' Person properties in a project
#'
#' Retrieve information about the item properties defined in the project (if any).
#' This will return a data.frame with one row per unique value of each item property. 
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data.frame with one row per unique value of each person property and a count of the nbr of persons it applies to.
#' The columns are: person_property, value, N
#'
get_person_properties <- function(db) {
  data = dbGetQuery(db, 'SELECT * FROM dxPersons;') 
  if(nrow(data) == 0 || ncol(data) ==1)
    return(data.frame())
  
  data %>% 
    gather(key='person_property', value='value', -.data$person_id) %>%
    group_by(.data$person_property, .data$value) %>%
    summarise(N=n())  %>%
	  ungroup()
}



#' Items in a project
#'
#' Retrieve all items that have been entered in the db
#' so far together with the item properties
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with column item_id and a column for each item property
#'
get_items <- function(db){
  dbGetQuery(db,'SELECT * FROM dxItems ORDER BY item_id;')
}


#' Persons in a project
#'
#' Retrieve all persons/respondents that have been entered in the db
#' so far together with their covariates
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with columns person_id and columns for each person_property/covariate
#'
get_persons <- function(db){
  dbGetQuery(db,'SELECT * FROM dxPersons ORDER BY person_id;')
}

#' Provide test scores
#'
#' Supplies the weighted sum of item scores for each person selected.
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @return A tibble with columns person_id, item_id, test_score
#' 
get_testscores<-function(dataSrc, predicate=NULL) {
  qtpredicate = eval(substitute(quote(predicate)))
  get_resp_data(dataSrc, qtpredicate, env=caller_env(), summarised=TRUE)$x %>%
    select(.data$person_id, .data$booklet_id, test_score=.data$sumScore) 
}


#' Test design
#'
#' Retrieve all items that have been entered in the db
#' so far by booklet and position in the booklet
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param format return format, see below
#' @param rows variable that defines the rows, ignored if format='long'
#' @param columns variable that defines the columns, ignored if format='long'
#' @param fill If set, missing values will be replaced with this value
#' @return A data.frame with the design. The contents depend on the rows, columns and format parameters
#'  if `format` is `'long'` a data.frame with columns: booklet_id, item_id, item_position
#'  if `format` is `'wide'` a data.frame with the rows defined by the `rows` parameter and 
#'  the columns by the `columns` parameter, with the remaining variable (i.e. item_id, booklet_id or item_position)
#'  making up the cells
#'
get_design <- function(db, 
                       format = c('wide','long'), 
                       rows = c('booklet_id','item_id','item_position'), 
                       columns = c('item_id','booklet_id','item_position'),
                       fill=NA)
{
  
  design = dbGetQuery(db,'SELECT booklet_id, item_id, item_position 
                            FROM dxBooklet_design
                              ORDER BY booklet_id, item_position;')
  format = match.arg(format)
  
  if(format == 'wide')
  {
    rows = match.arg(rows)
    columns = match.arg(columns)
    if(rows == columns) stop('rows may not be equal to columns')
    val_col = setdiff(c('booklet_id','item_id','item_position'), c(rows, columns))
    return(spread_(design, key_col = columns, value_col = val_col, fill = fill) %>% arrange_at(rows))
  } else
  {
    return(design)
  }

}




#' Test design as network
#'
#' Export the test design as an incidence matrix 
#' and a weight matrix
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
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
#' \donttest{
#' \dontrun{
#' dsgn = design_as_network(db)
#' # Check if design is connected
#' design_is_connected(dsgn)
#' }}
#'
design_as_network <- function(dataSrc, predicate = NULL, weights=c("items","responses")){
  
  # to do: the else catchall works but can be done more efficient in some cases. 
  # it would also be nice to support a design data frame instead of a responses dataframe in some way
  # mayb we can extend things yet further to do something with balanced designs and a priori designs
  # Left all of this for a future update
  
  w = match.arg(weights)
  qtpredicate = eval(substitute(quote(predicate)))
  
  if(is.null(qtpredicate) & inherits(dataSrc,'DBIConnection'))
  {
    if(w == 'items')
    {
      design = dbGetQuery(dataSrc, 'SELECT booklet_id, item_id FROM dxBooklet_design;')
    } else
    {
      # not entirely necessary to do a left join because booklet design can currently only be entered with respons data
      # but who knows we might want to separate that in the future to test for connectedness a priori
      design = dbGetQuery(dataSrc, 
                          'WITH pcount AS (SELECT booklet_id, item_id, COUNT(*) AS n FROM dxResponses GROUP BY booklet_id, item_id)
                          SELECT booklet_id, item_id, COALESCE(n,0) AS n_persons
                            FROM dxBooklet_design LEFT OUTER JOIN pcount USING(booklet_id, item_id);')
    }
  } else
  {
    if(w == 'items')
    {
      design = get_resp_data(dataSrc, qtpredicate, env=caller_env())$design
    } else
    {
      design = get_resp_data(dataSrc, qtpredicate, env=caller_env())$x  %>%
        group_by(.data$booklet_id, .data$item_id) %>%
        summarise(n_persons = n()) %>%
        ungroup()
    } 
  }
  
  if(length(unique(design$booklet_id)) < 2) stop("This makes sense only if you have at least two booklets")
  
  im = as.matrix(table(design$item_id, design$booklet_id))
  wm = crossprod(im, im)
  diag(wm) = 0
  if (w=="responses") {
    b = design %>% group_by(.data$booklet_id) %>% slice(1) %>% ungroup()
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
#' \donttest{
#' \dontrun{
#' # as an example, turn off some your booklets and see if you are
#' # still left with a connected design
#' 
#' dsgn = design_as_network(db, !(booklet_id %in% c('b1','b3','b4')))
#' design_is_connected(dsgn)
#' }
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




