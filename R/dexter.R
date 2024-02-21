utils::globalVariables(".")




#' Dexter: data analyses for educational and psychological tests.
#' 
#' Dexter provides a comprehensive solution for managing and analyzing educational test data.
#' 
#' The main features are:
#' 
#' \itemize{
#' \item project databases providing a structure for storing data about persons, items, responses and booklets.
#' \item methods to assess data quality using Classical test theory and plots.
#' \item CML calibration of the extended nominal response model and interaction model.
#' }
#' 
#' To learn more about dexter, start with the vignettes: `browseVignettes(package="dexter")`  
#' 
#' Dexter uses the following global options
#' \itemize{
#' \item `dexter.use_tibble` return tibbles instead of data.frames, defaults to FALSE
#' \item `dexter.progress` show progress bars, defaults to TRUE in interactive sessions
#' \item `dexter.max_cores` set a maximum number of cores that dexter will use, defaults to the minimum of `Sys.getenv("OMP_THREAD_LIMIT")` and
#' `getOption("Ncpus")`, otherwise unlimited.
#' }
#' 
"_PACKAGE"



#' Start a new project
#'
#' Imports a complete set of scoring rules and starts a new project (database)
#'
#'
#'
#' @param rules A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
#' The order is not important but spelling is. Any other columns will be ignored.
#' @param db_name A string specifying a filename
#' for a new sqlite database to be created. If this name does not
#' contain a path, the file will be created in the work
#' directory. Any existing file with the same name will be overwritten. For an in-memory database
#' you can use the string \code{":memory:"}. A connection object is also allowed.
#' @param person_properties An optional list of person properties. Names should correspond to person_properties intended to be used in the project.
#' Values are used as default (missing) values. The datatype will also be inferred from the values.
#' Known person_properties will be automatically imported when adding response data with \code{\link{add_booklet}}. 
#' @return a database connection object.
#' @details This package only works with closed items (e.g. likert, MC or possibly short answer)
#' it does not score any open items.
#' The first step to creating a project is to import an exhaustive list of all items and
#' all admissible responses, along with the score that any of the latter will be given.
#' Responses may be integers or strings but they will always be treated as strings.
#' Scores must be integers, and the minimum score for an item must be 0.
#' When inputting data, all responses not specified in the rules can optionally be treated as
#' missing and ultimately scored 0, but it is good style to include the missing
#' responses in the list. NA values will be treated as the string "NA"'.
#'
#' @examples
#'\donttest{
#' head(verbAggrRules)
#' db_name = tempfile(fileext='.db')
#' db = start_new_project(verbAggrRules, db_name, 
#'                        person_properties = list(gender = "unknown"))
#' }
#' 
start_new_project = function(rules, db_name="dexter.db", person_properties = NULL) 
{
  check_df(rules, c("item_id", "response", "item_score"))
  check_list(person_properties, nullable=TRUE)
  
  if(!is.null(person_properties))
      dbCheck_reserved_colnames(names(person_properties))

  rules = rules[, c("item_id", "response", "item_score")]
  rules$response = as.character(rules$response)
  rules$response[is.na(rules$response)] = 'NA'
  rules$item_id = as.character(rules$item_id)

  if(any(as.numeric(rules$item_score) %% 1 != 0))
    stop('only integer scores are allowed')
  
  rules$item_score = as.integer(rules$item_score)
  
  if(any(is.na(rules$item_id)) || any(is.na(rules$item_score)))
    stop("The item_id and item_score columns may not contain NA values")
  
  if(NROW(rules)>0)
  {
  	sanity = rules |>
  		group_by(.data$item_id) |>
  		summarise(less_than_two_scores = n_distinct(.data$item_score)<2,
  				duplicated_responses = any(duplicated(.data$response)),
  				min_score_not_zero = min(.data$item_score)>0) |>
  		filter(.data$less_than_two_scores | .data$duplicated_responses | .data$min_score_not_zero)

    if (nrow(sanity)>0)
  	{
        message("There were problems with your scoring rules.\nCheck the output below for possible reasons.\n")
        print(as.data.frame(sanity))
        stop('Scoring rules are not valid')
  	}
  } 

  if (inherits(db_name,'character'))
  {
      check_string(db_name)
      if (file.exists(db_name) && !suppressWarnings({file.remove(db_name)}))
        stop('Could not overwrite file: ', db_name)
      
      db = dbConnect(SQLite(), db_name)
  } else if(!is_db(db_name))
  {
      stop('argument db_name must be a filename or a DBIconnection')
  }
  dbTransaction(db,
  {
      project_CreateTables(db, person_properties)
      
      dbExecute_param(db,'INSERT INTO dxitems(item_id) VALUES(:id);', 
                      tibble(id=unique(rules$item_id)))
      dbExecute_param(db,
                'INSERT INTO dxscoring_rules(item_id, response, item_score) 
                          VALUES(:item_id, :response, :item_score);', 
				        rules)
  },
  on_error = function(e){dbDisconnect(db);stop(e)})

  db
}



#' Open an existing project
#'
#' Opens a database created by function \code{start_new_project}
#'
#'
#' @param db_name The name of the database to be opened.
#' @return a database connection object
#'
open_project = function(db_name="dexter.db") 
{
  check_string(db_name)
  if (!file.exists(db_name)) 
    stop("There is no such file")
    
  db = dbConnect(SQLite(), db_name)
  if(!dbExistsTable(db,'dxitems'))
  {
    dbDisconnect(db)
    stop('Sorry, this does not appear to be a Dexter database.')
  }
  dbExecute(db,'pragma foreign_keys=1;')
  return(db)
}

#' Close a project
#'
#' This is just an alias for \code{DBI::dbDisconnect(db)}, included for completeness
#'
#' @param db connection to a dexter database
#'
close_project = function(db) dbDisconnect(db) 


#' Derive scoring rules from keys
#'
#' For multiple choice items that will be scored as 0/1, derive the
#' scoring rules from the keys to the correct responses
#'
#'
#' @param keys  A data frame containing columns \code{item_id}, \code{noptions}, and
#' \code{key} See details.
#' @param include_NA_rule whether to add an option 'NA' (which is scored 0) to each item
#' @return A data frame that can be used as input to \code{start_new_project}
#' @details
#' This function might be useful in setting up the scoring rules when all items
#' are multiple-choice and scored as 0/1. 
#'
#' The input data frame must contain the exact id of each item, the number
#' of options, and the key. If the keys are all integers, it will be assumed that
#' responses are coded as 1 through noptions. If they are all  letters,
#' it is assumed that responses are coded as A,B,C,... All other cases result
#' in an error.
#'
keys_to_rules = function(keys, include_NA_rule = FALSE) 
{
  colnames(keys) = tolower(colnames(keys))
  check_df(keys, c('item_id','noptions','key'))
  keys = mutate_if(keys, is.factor, as.character)
  lwr = FALSE
  
  if (is.numeric(keys$key))
  {
    ABC=FALSE
  } 
  else 
  {
    if(all(keys$key %in% LETTERS))
    {
      ABC = TRUE
    } else if(all(keys$key %in% letters))
    {
      ABC=TRUE
      lwr=TRUE
      keys$key = toupper(keys$key)
    } else stop("You have inadmissible keys")
  }
  if (ABC) {
    m = match(keys$key, LETTERS)
    if (any(m>keys$noptions)) stop("You have out-of-range keys")
    r = keys |> 
      group_by(.data$item_id) |> 
      do({
      y = tibble(response=LETTERS[1:.$noptions[1]], item_score=0)
      y$item_score[match(.$key[1],LETTERS)] = 1
      y
    })
  } else {
    if (any(keys$key>keys$noptions)) stop("You have out-of-range keys")
    r = keys |> 
      group_by(.data$item_id) |> 
      do({
      y = tibble(response=as.character(1:.$noptions[1]), item_score=0)
      y$item_score[.$key[1]] = 1
      y
    })
  }
  r = ungroup(r)
  if(lwr)
  {
    r$response = tolower(r$response)
  }
  
  if(include_NA_rule){
    r = bind_rows(r, tibble(item_id=unique(r$item_id), response='NA', item_score=0))
  }
  r$item_score = as.integer(r$item_score)
  df_format(r)
}



#' Add or modify scoring rules
#' 
#' It is occasionally necessary to alter or add a scoring rule, e.g. in case of a key error. 
#' This function offers the possibility to do so and also allows you to add new items to your project
#' 
#' @param db a connection to a dexter project database
#' @param rules A data frame with columns \code{item_id}, \code{response}, and \code{item_score}.
#' The order is not important but spelling is. Any other columns will be ignored. See details
#' @return If the scoring rules pass a sanity check, a small summary of changes is printed and nothing is returned.
#' Otherwise this function returns a data frame listing the problems found, with 4 columns:
#' \describe{
#' \item{item_id}{id of the problematic item}
#' \item{less_than_two_scores}{if TRUE, the item has only one distinct score}
#' \item{duplicated_responses}{if TRUE, the item contains two or more identical response categories}
#' \item{min_score_not_zero}{if TRUE, the minimum score of the item was not 0}}
#' @details 
#' The rules should contain all rules that you want to change or add. This means that in case of a key error
#' in a single multiple choice question, you typically have to change two rules.
#' @examples 
#'\dontrun{\donttest{
#' # given that in your dexter project there is an mc item with id 'itm_01', 
#' # which currently has key 'A' but you want to change it to 'C'.
#' 
#' new_rules = data.frame(item_id='itm_01', response=c('A','C'), item_score=c(0,1))
#' touch_rules(db, new_rules)
#' }}
#' 
touch_rules = function(db, rules)
{
  check_df(rules, c('item_id','response','item_score'))
  check_db(db)
  
  if(any(is.na(rules$item_id)) || any(is.na(rules$item_score)))
    stop("The item_id and item_score columns may not contain NA values")
  
  rules = rules[, c("item_id", "response", "item_score")]
  rules$response = as.character(rules$response)
  rules$response[is.na(rules$response)] = 'NA'
  rules$item_id = as.character(rules$item_id)
  
  if(any(rules$item_score %% 1 != 0))
    stop('only integer scores are allowed')
  
  rules$item_score = as.integer(rules$item_score)
  
  # check if same options are supplied multiple times, this is the only check not done in conjunction
  # with the existing rules so it has to be gotten out of the way
  items_with_duplicate_options = unique(rules[duplicated(select(rules,'item_id', 'response')),]$item_id)

  existing_rules = dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxscoring_rules;')
  # remove the no-ops
  rules = dplyr::setdiff(rules, existing_rules)
  
  existing_opts = existing_rules |> select(-'item_score')
  
  # new items or responses
  new_rules = rules |> anti_join(existing_opts, by=c('item_id','response'))
  
  # existing items or responses but with new scores
  amended_rules = rules |> inner_join(existing_opts, by=c('item_id','response'))
  
  # to judge the validity of the new rules, we have to look at them in combination
  # with the rules in the db that will not be changed
  sanity = rules |> 
    dplyr::union(existing_rules |> 
                   anti_join(rules, by=c('item_id','response'))
                 ) |>
    group_by(.data$item_id) |>
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
      new_items = setdiff(new_rules$item_id, dbGetQuery(db, 'SELECT item_id FROM dxitems;')$item_id)
      if(length(new_items)>0) 
        dbExecute_param(db,'INSERT INTO dxitems(item_id) VALUES(:id);', tibble(id=new_items))
      
      dbExecute_param(db,'INSERT INTO dxscoring_rules(item_id, response, item_score) 
							                            VALUES(:item_id, :response, :item_score);', 
				select(new_rules, 'item_id', 'response', 'item_score'))
    }
    if(nrow(amended_rules)>0) 
    {
      dbExecute_param(db,'UPDATE dxscoring_rules SET item_score=:item_score 
                      WHERE item_id=:item_id AND response=:response;', 
                select(amended_rules, 'item_id', 'response', 'item_score'))
    }
  })
  cat(paste0('\nrules_changed: ', nrow(amended_rules), '\nrules_added: ', nrow(new_rules)),'\n')
}


#' Add response data to a project
#'
#' Add item response data in long or wide format.
#'
#'
#' @param db a connection to a dexter database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param x A data frame containing the responses and, optionally,
#' person_properties. The data.frame should have one row per respondent and the column names should 
#' correspond to the item_id's in the rules or the names of the person_properties. See details.
#' @param booklet_id A (short) string identifying the test form (booklet)
#' @param data response data in normalized (long) format. Must contain columns \code{person_id}, \code{booklet_id},
#'  \code{item_id} and \code{response} and optionally \code{item_position} 
#' (useful if your data contains new booklets, see details)
#' @param design data.frame with columns booklet_id, item_id and optionally item_position specifying the design of any 
#' _new_ booklets in your data.
#' @param missing_value value to use for responses in missing rows in your data, see details
#' @param auto_add_unknown_rules  If FALSE (the default), an error will be generated if 
#' one or more responses do not appear in the scoring rules. If TRUE, unknown responses will be
#' assumed to have a score of 0 and will be added to your scoring rules
#' 
#' @return A list with information about the recent import.
#' 
#' @details It is a common practice to keep response data in tables where each row 
#' contains the responses from a single person. \code{add_booklet} is provided to input
#' data in that form, one booklet at a time. 
#'
#' If the dataframe \code{x} contains a variable named \code{person_id} this variable 
#' will be used to identify unique persons. It is assumed that a single person will only 
#' make a single booklet once, otherwise an error will be generated. 
#' 
#' If a person_id is not supplied, dexter will generate unique person_id's for each row of data.  
#'
#' Any column whose name has an exact match in the scoring rules inputted with
#' function \code{start_new_project} will be treated as an item; any column whose name has an 
#' exact match in the person_properties will be treated as a person property. If a name matches both
#' a person_property and an item_id, the item takes precedence. Columns other than items, person properties 
#' and person_id will be ignored.
#' 
#' 
#' \code{add_response_data} can be used to add data that is already normalized. This function takes a 
#' data.frame in long format with columns \code{person_id}, \code{booklet_id}, 
#' \code{item_id} and \code{response} such as can usually be found in databases for example. 
#' For booklets that are not already known in your project, you need to specify the design via the \code{design} argument.
#' Failure to do so will result in an error. Responses to items that should be there according to the design but which do not have a corresponding
#' row in \code{data} will be added with \code{missing_value} used for the response. If this missing value is not defined in your scoring rules 
#' and \code{auto_add_unknown_rules} is set to FALSE, this will lead to an error message.
#' 
#' 
#' Note that responses are always treated as strings (in both functions), and \code{NA}
#' values are transformed to the string \code{"NA"}.
#' 
#' @examples 
#' db = start_new_project(verbAggrRules, ":memory:", 
#'                        person_properties=list(gender="unknown"))
#' head(verbAggrData)
#' add_booklet(db, verbAggrData, "agg")      
#' 
#' close_project(db)
#' 
add_booklet = function(db, x, booklet_id, auto_add_unknown_rules = FALSE) {
  
  check_df(x)
  check_db(db)
  x = mutate_if(x, is.factor, as.character) 
  
  person_properties = intersect(dbListFields(db, 'dxpersons'), tolower(names(x)))
  person_properties = person_properties[person_properties != 'person_id']
  
  design = tibble(booklet_id = booklet_id, item_id = names(x), col_order=c(1:ncol(x))) |>
    inner_join(dbGetQuery(db, "SELECT item_id FROM dxitems;"), by='item_id') |>
    mutate(item_position = dense_rank(.data$col_order)) |>
    select(-'col_order')

  if(nrow(design) == 0) stop('None of the column names in x correspond to known items in your project')

  dbTransaction(db,{
    out = list()
    is_existing_booklet = dbExists(db,'SELECT 1 FROM dxbooklets WHERE booklet_id=:b;',tibble(b=booklet_id))
    if (is_existing_booklet) 
    {
      if(!df_identical(dbGetQuery_param(db,
              'SELECT item_id FROM dxbooklet_design WHERE booklet_id=:b ORDER BY item_position;', 
              tibble(b=booklet_id)),
           tibble(item_id = design$item_id)))
      {
        stop_("There is already a booklet with this ID which has different items or a different item order")
      }
    } else
    {  
      dbExecute_param(db,'INSERT INTO dxbooklets(booklet_id) VALUES(:booklet_id);',list(booklet_id=booklet_id))
      if('testpart_nbr' %in% dbListFields(db, 'dxbooklet_design'))
      {
        dbExecute_param(db,'INSERT INTO dxtestparts(booklet_id) VALUES(:booklet_id);',list(booklet_id=booklet_id))
      }
      dbExecute_param(db,'INSERT INTO dxbooklet_design(booklet_id, item_id, item_position) 
                          VALUES(:booklet_id,:item_id,:item_position);', 
				select(design, 'booklet_id', 'item_id', 'item_position'))
    }
                  
    x$booklet_id = booklet_id
    if(!'person_id' %in% names(x))
    {
      x$person_id = dbUniquePersonIds(db,nrow(x))
      new_people = x$person_id
      message("no column `person_id` provided, automatically generating unique person id's")
    } else
    {
      x$person_id = as.character(x$person_id)
      known_people = dbGetQuery(db,'SELECT person_id FROM dxpersons;')$person_id
      new_people = setdiff(x$person_id,known_people)
      if(anyDuplicated(x$person_id))
        stop_("The column `person_id` has duplicate values, this is not allowed.")

      if(length(new_people) < nrow(x) && is_existing_booklet)
      {
        known_admin = dbGetQuery_param(db,"SELECT person_id FROM dxadministrations WHERE booklet_id=:b;",
                                       tibble(b=booklet_id))$person_id
        if(length(intersect(known_admin, x$person_id))>0)
          stop_("One or more person_id's in your data overlap with persons already imported who made the same booklet. This is not allowed. ")
      }
    }
                  
    dbExecute_param(db,'INSERT INTO dxpersons(person_id) VALUES(:person_id);', tibble(person_id=new_people))
    
    dbExecute_param(db,'INSERT INTO dxadministrations(person_id,booklet_id) VALUES(:person_id,:booklet_id);', 
              select(x, 'person_id', 'booklet_id'))
    
   
    responses = x |>
      select(all_of(c(design$item_id, "booklet_id", "person_id"))) |>
      pivot_longer(all_of(design$item_id), values_drop_na = FALSE,
                   names_to='item_id', values_to='response',values_transform=as.character)
    
    responses$response[is.na(responses$response)] = 'NA'
                  
    new_rules = anti_join(responses, dbGetQuery(db, "SELECT item_id, response FROM dxscoring_rules;"), 
                          by=c('item_id','response')) |>
      distinct(.data$item_id, .data$response)
    
    if (nrow(new_rules)>0 && auto_add_unknown_rules) 
    {
		  dbExecute_param(db,
				'INSERT INTO dxscoring_rules(item_id,response,item_score) VALUES(:item_id,:response,0);',
				select(new_rules, 'item_id', 'response'))
      out$zero_rules_added = new_rules
    } else if(nrow(new_rules)>0) 
    {
      message('The following responses are not in your rules (showing first 30):\n')
      print(head(as.data.frame(new_rules), 30), row.names=FALSE)
      stop_('unknown responses')
    }
    dbExecute_param(db,'INSERT INTO dxresponses(booklet_id,person_id,item_id,response) 
                                VALUES(:booklet_id,:person_id,:item_id,:response);', 
					select(responses, 'booklet_id', 'person_id', 'item_id', 'response'))
            
    # make this report before we mutilate the colnames  
    columns_ignored = setdiff(names(x), c(design$item_id,'person_id','item_id','booklet_id') )
    columns_ignored = columns_ignored[!tolower(columns_ignored) %in% person_properties]
                  
    # add person_properties
    if(length(person_properties)>0)
    {
      colnames(x) = tolower(colnames(x))
      dbExecute_param(db,paste0('UPDATE dxpersons SET ',paste0(person_properties,'=:',person_properties,collapse=','),' WHERE person_id=:person_id;'),
                               x[,c(person_properties,'person_id')])
    }
  }) 
  
  out$items = design$item_id
  out$person_properties = person_properties
  if( length(columns_ignored) > 0)
    out$columns_ignored = columns_ignored

  out
}

#' @rdname add_booklet 
add_response_data = function(db, data, design=NULL, missing_value = 'NA', auto_add_unknown_rules = FALSE)
{
  colnames(data) = tolower(colnames(data))
  
  check_df(data, c('item_id', 'person_id', 'response','booklet_id'))
  check_db(db)
  
  data = data[,c('item_id', 'person_id', 'response','booklet_id')]
  data = mutate_all(ungroup(data), as.character)
  
  missing_value = as.character(missing_value)
  check_string(missing_value)
  
  db_booklets = dbGetQuery(db, "SELECT booklet_id FROM dxbooklets;")[,1]
  
  if(is.null(design))
  {
    new_bk = setdiff(data$booklet_id, db_booklets)
    if(length(new_bk) > 0)
    {
      message('Unknown booklets:')
      print(new_bk)
      stop_('`data` contains unknown booklets. You have to specify any new booklets via the design argument.')
    }
  } else
  {
    check_df(design,c('booklet_id','item_id'))
    design = ungroup(design)
    design$booklet_id = as.character(design$booklet_id)
    design$item_id = as.character(design$item_id)
    if('item_position' %in% colnames(design))
    {
      design$item_position = as.integer(design$item_position)
    } else
    {
      design = design |>
        group_by(.data$booklet_id) |>
        mutate(item_position = row_number()) |>
        ungroup()
    }
    design = design[,c('booklet_id','item_id','item_position')]
    
    if(nrow(design) > n_distinct(design$booklet_id, design$item_id))
      stop_('booklet_id, item_id combination in `design` is not unique')
    if(nrow(design) > n_distinct(design$booklet_id, design$item_position))
      stop_('booklet_id, item_position combination in `design` is not unique')
    
    unknown_items = setdiff(design$item_id,dbGetQuery(db,'SELECT item_id FROM dxitems;')[,1])
    if(length(unknown_items)>0)
    {
      message('Unknown items:')
      print(unknown_items)
      stop_('`design` contains unknown items. You have to specify any new items and scoring rules using the function touch_rules.')
    }
    
    
    if(any(design$booklet_id %in% db_booklets))
    {
      check_dsg_db = get_design(db) |>
        semi_join(design,by='booklet_id') |>
        arrange(.data$booklet_id, .data$item_id)
      
      check_dsg = design |>
        filter(.data$booklet_id %in% db_booklets) |>
        arrange(.data$booklet_id, .data$item_id)
      
      
      # ignore item position in check
      if(nrow(check_dsg) == nrow(check_dsg_db) && 
         all(check_dsg_db$booklet_id == check_dsg$booklet_id) && 
         all(check_dsg_db$item_id == check_dsg$item_id))
      {
        design = anti_join(design, check_dsg, by='booklet_id')
        if(NROW(design) == 0) design = NULL
      } else
      {
        stop_("Design contains booklets that are already defined in the database but don't have the same items.")
      }
    }
  }
  if(n_distinct(data$person_id,data$booklet_id,data$item_id) < nrow(data))
    stop_("The combination of the columns person_id, booklet_id, item_id in `data` is not unique")
  
  
  msg = NULL
  
  dbTransaction(db,{ 
    # update design
    if(!is.null(design))
    {
      dbExecute_param(db,'INSERT INTO dxbooklets(booklet_id) VALUES(:booklet_id);',distinct(design, .data$booklet_id))
      dbExecute_param(db,'INSERT INTO dxbooklet_design(booklet_id, item_id, item_position) VALUES(:booklet_id,:item_id,:item_position);',design)
    }
    
    # check data
    data_design = distinct(data, .data$booklet_id, .data$item_id)
    design = semi_join(get_design(db), data_design, by='booklet_id') |> 
      select('booklet_id', 'item_id')
    
    invalid_bk_item = anti_join(data_design, design, by=c('booklet_id','item_id'))
    if(NROW(invalid_bk_item) > 0)
    {
      message('Unknown booklet, item combinations (showing first 10)')
      print(head(invalid_bk_item,10))
      stop_('Your data contains booklet,item combinations that should not occur according to the design.')
    }
    
    # fill in the missings
    administrations = distinct(data,.data$booklet_id, .data$person_id)
    
    data = administrations |>
      inner_join(design, by='booklet_id', relationship = "many-to-many") |>
      left_join(data, by=c('person_id','booklet_id','item_id'))
    
    NA_cnt = sum(is.na(data$response))
    if(NA_cnt > 0)
    {
      data$response[is.na(data$response)] = missing_value
      msg = c(msg,sprintf('%i missing responses replaced by "%s"', NA_cnt, missing_value))
    }
    
    # check rules
    new_rules = distinct(data, .data$item_id, .data$response) |>
      anti_join(get_rules(db), by=c('item_id','response'))
    
    if(NROW(new_rules)>0)
    {
      unknown_items = distinct(new_rules, .data$item_id) |>
        anti_join(dbGetQuery(db, 'SELECT item_id FROM dxitems;'), by='item_id')
      
      if(NROW(unknown_items) > 0)
      {
        message('The following items are not known in your project:')
        print(unknown_items$item_id) 
        stop("encountered item_id's not defined in your project")
      }
      
      if(auto_add_unknown_rules)
      {
        dbExecute_param(db,"INSERT INTO dxscoring_rules(item_id, response, item_score) VALUES(:item_id, :response, 0);",new_rules)
        msg = c(msg,sprintf('%i scoring rules with 0 score added to the rules', nrow(new_rules)))
      } else
      {
        if(all(new_rules$response == missing_value))
        {
          message(sprintf('%i missing responses replaced by "%s", but this value is not known in the scoring rules.\nType ?add_response_data for help.',
                          NA_cnt, missing_value))
        } else
        {
          message('Unknown responses (showing first 10):')
          print(head(new_rules,10))
        }
        stop_("Unknown responses")
      }
    }
    
    # everything seems alright
    # add persons
    persons = distinct(administrations, .data$person_id)
    new_persons = anti_join(persons, dbGetQuery(db, 'SELECT person_id FROM dxpersons;'), by='person_id')
    if(nrow(new_persons) > 0)
      dbExecute_param(db, 'INSERT INTO dxpersons(person_id) VALUES(:person_id);', new_persons)
    
    # add administrations
    existing_admin = inner_join(administrations, 
                                dbGetQuery(db,'SELECT person_id, booklet_id FROM dxadministrations;'), 
                                by=c('person_id','booklet_id'))
    
    if(nrow(existing_admin) > 0)
    {
      message('The following person-booklet combination have already been entered into the project (showing first 20)')
      head(existing_admin, 20) |> 
        as.data.frame() |> 
        print(row.names=FALSE)
      stop('double administrations')
    }
    
    dbExecute_param(db, 'INSERT INTO dxadministrations(person_id, booklet_id) VALUES(:person_id, :booklet_id);', 
                    administrations)
    
    # add data
    n = dbExecute_param(db, 
                        'INSERT INTO dxresponses (person_id, booklet_id, item_id, response) VALUES(:person_id, :booklet_id, :item_id, :response);', 
                        data)
  })
  msg = c(msg,paste(n,'responses imported.\n'))
  cat(paste(msg,collapse='\n'))
}



#' Add item properties to a project
#'
#' Add, change or define item properties in a dexter project
#'
#'
#' @param db a connection to a dexter database, e.g. the output of \code{start_new_project}
#' or \code{open_project}
#' @param item_properties A data frame containing a column item_id (matching item_id's already defined in the project)
#'  and 1 or more other columns with item properties (e.g. item_type, subject)
#' @param default_values a list where the names are item_properties and the values are defaults.
#' The defaults will be used wherever the item property is unknown.
#' @return nothing
#' @details When entering response data in the form of a rectangular person x item
#' table, it is easy to provide person properties but practically impossible
#' to provide item properties. This function provides a possibility to do so.
#' 
#' Note that is is not possible to add new items with this function, 
#' use \code{\link{touch_rules}} if you want to add new items to your project.
#' 
#' @seealso \code{\link{fit_domains}}, \code{\link{profile_plot}} for
#'  possible uses of item_properties
#'
#' @examples 
#' \dontrun{\donttest{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' head(verbAggrProperties)
#' add_item_properties(db, verbAggrProperties)
#' get_items(db) 
#' 
#' close_project(db)
#' }}
#'
add_item_properties = function(db, item_properties=NULL, default_values=NULL) {
  
  check_db(db)
  check_df(item_properties, 'item_id', nullable=TRUE)
  check_list(default_values, nullable=TRUE)
  
  # for convenience we include the item_id as a property
  existing_item_properties = dbListFields(db, 'dxitems') 
  
  dbTransaction(db, 
  {
    if(!is.null(default_values))
    {
      stopifnot(is.list(default_values))
      names(default_values) = dbValid_colnames(names(default_values))
      
      new_prop_names = setdiff(names(default_values), existing_item_properties)
      dbCheck_reserved_colnames(new_prop_names)
      
      for(prop_name in new_prop_names)
      {
        dbExecute(db, 
                  paste0("ALTER TABLE dxitems ADD COLUMN ",
                         prop_name, 
                         sql_col_def(default_values[[prop_name]], TRUE, db),';'))
      }
      cat(paste(length(new_prop_names),'new item_properties defined\n'))
      existing_item_properties = c(existing_item_properties, new_prop_names)
    }
    if(!is.null(item_properties))
    {
      colnames(item_properties) = dbValid_colnames(colnames(item_properties))
      
      item_properties = item_properties |>
        mutate(item_id = as.character(.data$item_id)) |>
        mutate_if(is.factor, as.character) 
      
      if(inherits(db,'SQLiteConnection'))
      {
        item_properties = item_properties |>
          mutate_if(is.date, format, "%Y-%m-%d") |>
          mutate_if(is.time, format, "%Y-%m-%d %H:%M:%S")
      }

      if(!('item_id' %in% colnames(item_properties)))
        stop('item_properties needs to have a column item_id')
      
      new_prop_names = setdiff(colnames(item_properties), existing_item_properties)
      dbCheck_reserved_colnames(new_prop_names)
      
      for(prop_name in new_prop_names)
      {
        dbExecute(db, 
                paste0("ALTER TABLE dxitems ADD COLUMN ", 
                       prop_name, 
                       sql_data_type(pull(item_properties, prop_name)),";"))
      }
      pnames = names(item_properties)[names(item_properties)!='item_id']
      
      n = dbExecute_param(db,
            paste0('UPDATE dxitems SET ',paste0(pnames,'=:',pnames,collapse=', '),
                      ' WHERE item_id=:item_id;'),
            item_properties)
      cat(paste(length(pnames), 'item properties for', n, 'items added or updated\n'))
    }
  })
}



#' Add person properties to a project
#'
#' Add, change or define person properties in a dexter project. Person properties defined here will 
#' also be automatically imported with \code{\link{add_booklet}}
#'
#'
#' @param db a connection to a dexter database, e.g. the output of \code{start_new_project}
#' or \code{open_project}
#' @param person_properties A data frame containing a column person_id and 1 or more other columns with 
#' person properties (e.g. education_type, birthdate)
#' @param default_values a list where the names are person_properties and the values are defaults.
#' The defaults will be used wherever the person property is unknown.
#' 
#' @details
#' Due to limitations in the sqlite database backend that we use, the default values for a person property 
#' can only be defined once for each person_property
#' 
#' @return nothing
add_person_properties = function(db, person_properties = NULL, default_values = NULL)
{
  check_db(db)
  check_df(person_properties, 'person_id', nullable=TRUE)
  check_list(default_values, nullable=TRUE)
  
  dbTransaction(db, 
  { 
    existing_props = dbListFields(db, 'dxpersons')
    
    if(!is.null(default_values))
    {
      names(default_values) = dbValid_colnames(names(default_values))
      
      dbCheck_reserved_colnames(names(default_values))
  
      for(prop_name in setdiff(names(default_values), existing_props))
      {
        dbExecute(db, paste0("ALTER TABLE dxpersons ADD COLUMN ",prop_name, sql_col_def(default_values[[prop_name]], TRUE, db),';'))
      }
      cat(paste(length(setdiff(names(default_values), existing_props)),"new person_properties defined\n"))
      existing_props = union(existing_props, names(default_values))
    }
    if(!is.null(person_properties))
    {
      person_properties$person_id = as.character(person_properties$person_id)
      colnames(person_properties) = dbValid_colnames(colnames(person_properties))

      dbCheck_reserved_colnames(setdiff(colnames(person_properties), existing_props))
      
      if(inherits(db,'SQLiteConnection'))
      {
        person_properties = person_properties |>
          mutate_if(is.date, format, "%Y-%m-%d") |>
          mutate_if(is.time, format, "%Y-%m-%d %H:%M:%S")
      }
      
      
      lapply(setdiff(colnames(person_properties), existing_props), function(prop_name)
      {
        dbExecute(db, paste0("ALTER TABLE dxpersons ADD COLUMN ",
                             prop_name, 
                             sql_col_def(pull(person_properties, prop_name), FALSE, db),';'))
      })
      
      pnm = setdiff(colnames(person_properties),'person_id')
      
      n = dbExecute_param(db, 
                    paste('UPDATE dxpersons SET', paste0(pnm, '=:', pnm, collapse = ','),
                              'WHERE person_id=:person_id;'),
                    person_properties)
      
      cat(paste(ncol(person_properties) - 1, 'person properties for', n, 'persons added or updated\n'))
    }
  })
}




#' Get scoring rules
#' 
#' Retrieve the scoring rules currently present in the dexter project db
#' 
#' @param db a connection to a Dexter database
#' @return data.frame of scoring rules containing columns: item_id, response, item_score
#' 
get_rules = function(db)
{
  check_db(db)
  df_format(dbGetQuery(db, 'SELECT item_id, response, item_score FROM dxscoring_rules ORDER BY item_id, response;'))
}


#' Booklets entered in a project
#'
#' Retrieve information about the booklets entered in the db so far
#'
#' @param db a connection to a dexter database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with columns: booklet_id, n_persons, n_items and booklet_max_score. 
#' booklet_max_score gives the maximum theoretically possible score according to the scoring rules
#'
get_booklets = function(db) {
  check_db(db)
  bk = dbGetQuery(db, 'SELECT dxbooklets.*, n_items, n_persons FROM dxbooklets 
                        INNER JOIN dxbooklet_stats USING(booklet_id) ORDER BY booklet_id;')
  if(!'booklet_max_score' %in% colnames(bk))
  {
    ms = dbGetQuery(db,
      "WITH I AS (SELECT item_id, MAX(item_score) AS ms
        FROM dxscoring_rules GROUP BY item_id)
      
      SELECT booklet_id, SUM(ms) AS booklet_max_score
        FROM dxbooklet_design INNER JOIN I USING(item_id)
          GROUP BY booklet_id;")
    bk = inner_join(bk,ms,by='booklet_id')
  }
  df_format(bk)
}




#' Variables that are defined in the project
#' 
#' Inspect the variables defined in your dexter project and their datatypes
#' 
#' @param db a dexter project database
#' @return a data.frame with name and type of the variables defined in your dexter project
#' @details 
#' The variables in Dexter consist of the item properties and person properties you specified
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
  check_db(db)
  lapply(c('dxitems','dxbooklets','dxbooklet_design','dxscoring_rules',
           'dxpersons','dxadministrations','dxresponses'),
         function(tbl)
         {
           res = dbSendQuery(db,paste('SELECT * FROM',tbl,'WHERE 0=1;'))
           r = dbColumnInfo(res)
           dbClearResult(res)
           r$originates = substring(tbl,3)
           return(r)
         }) |> 
    bind_rows() |> 
    distinct(.data$name,.keep_all=TRUE) |>
    filter(.data$name != 'testpart_nbr') |>
    arrange(.data$originates)
}




#' Items in a project
#'
#' Retrieve all items that have been entered in the db
#' so far together with the item properties
#'
#'
#' @param db a connection to a dexter database, e.g. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with column item_id and a column for each item property
#'
get_items = function(db)
{
  if(is.list(db) && !is.null(db$inputs))
    return(df_format(tibble(item_id=as.character(db$inputs$ssI$item_id))))

  check_db(db)
  df_format(dbGetQuery(db,'SELECT * FROM dxitems ORDER BY item_id;'))
}


#' Persons in a project
#'
#' Retrieve all persons/respondents that have been entered in the db
#' so far together with their properties
#'
#'
#' @param db a connection to a dexter database, e.g. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame with columns person_id and columns for each person_property
#'
get_persons = function(db){
  check_db(db)
  df_format(dbGetQuery(db,'SELECT * FROM dxpersons ORDER BY person_id;'))
}

#' Get test scores
#'
#' Supplies the sum of item scores for each person selected.
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @return A tibble with columns person_id, item_id, booklet_score
#' 
get_testscores = function(dataSrc, predicate=NULL) 
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  get_resp_data(dataSrc, qtpredicate, summarised=TRUE, env=env)$x |>
    mutate_if(is.factor, as.character) |>
    df_format()
}


#' Test design
#'
#' Retrieve all items that have been entered in the db
#' so far by booklet and position in the booklet
#'
#'
#' @param dataSrc a dexter database or any object form which a design can be inferred
#' @param format return format, see below
#' @param rows variable that defines the rows, ignored if format='long'
#' @param columns variable that defines the columns, ignored if format='long'
#' @param fill If set, missing values will be replaced with this value, ignored if format='long'
#' @return A data.frame with the design. The contents depend on the rows, columns and format parameters
#'  if \code{format} is \code{'long'} a data.frame with columns: booklet_id, item_id, item_position (if available)
#'  if \code{format} is \code{'wide'} a data.frame with the rows defined by the \code{rows} parameter and 
#'  the columns by the \code{columns} parameter, with the remaining variable (i.e. item_id, booklet_id or item_position)
#'  making up the cells
#'
get_design = function(dataSrc, 
                       format = c('long','wide'), 
                       rows = c('booklet_id','item_id','item_position'), 
                       columns = c('item_id','booklet_id','item_position'),
                       fill=NA)
{
  
  design = 
    if(is_db(dataSrc))
    {
      dbGetQuery(dataSrc,'SELECT * 
                            FROM dxbooklet_design
                                ORDER BY booklet_id, item_position;')
    } else if(inherits(dataSrc, 'list') && !is.null(dataSrc$inputs$design))
    {
      mutate_if(dataSrc$inputs$design[,c('booklet_id','item_id')], is.factor, as.character)
    } else if(inherits(dataSrc, 'data.frame'))
    {
      colnames(dataSrc) = tolower(colnames(dataSrc))
      dataSrc[,intersect(c('booklet_id','item_position','item_id'),colnames(dataSrc))] |>
        distinct() |>
        arrange_all()
    } else
    {
      get_resp_data(dataSrc)$design |>
        distinct() |>
        arrange_all()
    }
  
  
  format = match.arg(format)
  
  if(format == 'wide')
  {
    if(!'item_position' %in% colnames(design))
      stop('Cannot make wide format becasue item position is missing in input')

    rows = match.arg(rows)
    columns = match.arg(columns)
    if(rows == columns) stop('rows may not be equal to columns')
    val_col = setdiff(c('booklet_id','item_id','item_position'), c(rows, columns))

    design |>
      select('item_id','item_position','booklet_id') |>
      pivot_wider(names_from=columns, values_from=val_col, values_fill=fill, names_sort=TRUE) |>
      arrange(all_of(rows)) |>
      df_format()
        
  } else
  {
    df_format(design)
  }
}



# to~do: add a unit test

#' Information about the design
#'
#' This function is useful to inspect incomplete designs
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' 
#' 
#' @return
#' a list with the following components
#' \describe{
#' \item{design}{a data.frame with columns booklet_id, item_id, item_position, n_persons} 
#' \item{connected_booklets}{a data.frame with columns booklet_id, group; 
#' booklets with the same `group` are connected to each other.} 
#' \item{connected}{TRUE/FALSE indicating whether the design is connected or not} 
#' \item{testlets}{a data.frame with columns item_id and testlet; items within the same testlet 
#' always occur together in a booklet} 
#' \item{adj_matrix}{list of two adjacency matrices: *weighted_by_items* and *weighted_by_persons*; These matrices can be 
#' useful in visually inspecting the design using a package like *igraph*}
#' }
#' 
design_info = function(dataSrc, predicate = NULL)
{
  out = list()
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()

  #check_dataSrc(dataSrc) # removed, prevents design being inputted

  # parms objects
  if(inherits(dataSrc,'list') && !is.null(dataSrc$inputs$design))
    dataSrc = dataSrc$inputs$design

  if(inherits(dataSrc, 'data.frame'))
    colnames(dataSrc) = tolower(colnames(dataSrc))
  
  if(inherits(dataSrc, 'data.frame') && !('person_id'%in% colnames(dataSrc)))
  {
    # perhaps a design was supplied
    out$design = dataSrc |>
      distinct(.data$booklet_id, .data$item_id, .keep_all=TRUE)
    
    out$design = out$design[,intersect(colnames(out$design), c('booklet_id','item_id','item_position'))]
    out$design$n_persons = NA_integer_
  } else
  {
    check_dataSrc(dataSrc)
    if(is_db(dataSrc))
    {
      out$design = db_get_design(dataSrc, qtpredicate=qtpredicate, env=env)
    } else
    {
      out$design = get_resp_data(dataSrc, qtpredicate, env=env)$x |>
        count(.data$booklet_id, .data$item_id, name='n_persons')
    }
  }

  out$design = out$design |>
    arrange(.data$booklet_id, .data$item_id)  |> 
    mutate_if(is.factor, as.character) |>
    df_format()
      
    
  if(nrow(out$design) == 0 )
    stop("design is empty: no data selected")
    
 
  #testlets
  out$testlets = out$design |>
    mutate(bnr = dense_rank(.data$booklet_id)) |>
    group_by(.data$item_id) |>
    summarize(testlet = paste(sort(.data$bnr), collapse = ' ')) |>
    ungroup() |>
    mutate(testlet = dense_rank(.data$testlet)) |>
    arrange(.data$testlet, .data$item_id) |>
    df_format()
  
  
  # to~do: I don't see why the diagonals should be zero, check in igraph if this matters
  # network
  out$adj_matrix = list()
  
  items = as.matrix(table(out$design$item_id, out$design$booklet_id))
  out$adj_matrix$weighted_by_items = crossprod(items, items)
  mode(out$adj_matrix$weighted_by_items) = 'integer'
  diag(out$adj_matrix$weighted_by_items) = 0L
  
  
  b = out$design |> 
    distinct(.data$booklet_id, .keep_all=TRUE)
  ww = outer(b$n_persons, b$n_persons, "+")
  out$adj_matrix$weighted_by_persons = out$adj_matrix$weighted_by_items * ww    

  
  
  out$connected_booklets = tibble(booklet_id = colnames(out$adj_matrix$weighted_by_items), 
                                  group = ds_connected_groups(out$adj_matrix$weighted_by_items)) |> 
    df_format()
  out$connected = (max(out$connected_booklets$group) == 1)
  
  
  out
}



# datasets



#' Verbal aggression data
#' 
#' A data set of self-reported verbal behaviour in different frustrating
#' situations (Vansteelandt, 2000). The dataset also contains participants reported gender and scores on the 'anger' questionnaire.
#' 
#' 
#' 
#' @name verbAggrData
#' @docType data
#' @format A data set with 316 rows and 26 columns.
#' @keywords datasets
NULL


#' Scoring rules for the verbal aggression data
#' 
#' A set of (trivial) scoring rules for the verbal 
#' aggression data set
#' 
#' 
#' @name verbAggrRules
#' @docType data
#' @format A data set with 72 rows and 3 columns (item_id, response, item_score).
#' @keywords datasets
NULL


#' Item properties in the verbal aggression data
#' 
#' A data set of item properties related to the verbal
#' aggression data
#' 
#' 
#' @name verbAggrProperties
#' @docType data
#' @format A data set with 24 rows and 5 columns.
#' @keywords datasets
NULL


#' Rated data
#' 
#' A data set with rated data. A number of student performances are rated twice on several
#' aspects by independent judges. The ratings are binary and have been summed following
#' the theory discussed by Maris and Bechger (2006, Handbook of Statistics). Data are a
#' small subset of data collected on the State Exam Dutch as a second language for Speaking. 
#' 
#' 
#' 
#' @name ratedData
#' @docType data
#' @format A data set with 75 rows and 15 columns.
#' @keywords datasets
NULL


#' Scoring rules for the rated data
#' 
#' A set of (trivial) scoring rules for the rated data set
#' 
#' 
#' @name ratedDataRules
#' @docType data
#' @format A data set with 42 rows and 3 columns (item_id, response, item_score).
#' @keywords datasets
NULL


#' Item properties in the rated data
#' 
#' A data set of item properties related to the rated data. These are the aspects: 
#' IH = content, WZ = word choice and phrasing, and WK = vocabulary.
#' 
#' 
#' @name ratedDataProperties
#' @docType data
#' @format A data set with 14 rows and 2 columns: item_id and aspect
#' @keywords datasets
NULL


