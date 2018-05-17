

#' Selecting data
#' 
#' Extract data from a dexter database 
#' 
#' @param dataSrc a dexter project database or data.frame
#' @param predicate an expression to select data on
#' @param columns the columns you wish to select, can include all columns in the project, see: \code{\link{get_variables}}
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
get_responses = function(dataSrc, predicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) env = caller_env() 
  as.data.frame(get_responses_(dataSrc, eval(substitute(quote(predicate))), columns=columns, env=env)) 
}

get_responses_ = function(dataSrc, qtpredicate=NULL, columns=c('person_id','item_id','item_score'), env=NULL)
{
  if(is.null(env)) env=caller_env() 
  
  if(inherits(dataSrc,'data.frame'))
  {
    if(is.null(qtpredicate))
    {
      respData = dataSrc[,columns]
    } else
    {
      respData = dataSrc[eval_tidy(qtpredicate, data=dataSrc, env=env), columns]
    }
  } else
  {
    db = dataSrc
    columns = dbValid_colnames(columns)
    if(is.null(qtpredicate))
    {
      where = ''
      used_columns = columns
    } else
    {
      suppressWarnings({where = qtpredicate2where(qtpredicate, db, env)})
      used_columns = union(columns, tolower(all.vars(qtpredicate)))                      
    }
    # decide which tables we need to join since joining costs time
    used_columns = setdiff(used_columns, 
                           c(dbListFields(db,'dxResponses'), dbListFields(db,'dxScoring_rules')))
    cte = " dxResponses INNER JOIN dxScoring_rules USING(item_id, response)"
    if(length(intersect(dbListFields(db,'dxPersons'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxPersons USING(person_id)')
    if(length(intersect(dbListFields(db,'dxItems'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxItems USING(item_id)')
    if(length(intersect(dbListFields(db,'dxBooklets'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxBooklets USING(booklet_id)')
    if(length(intersect(dbListFields(db,'dxBooklet_design'),used_columns))>0) 
      cte = c(cte, 'INNER JOIN dxBooklet_design USING(booklet_id, item_id)')
    
    # can have unexpected behavior if person makes more than one testform
    respData = NULL
    tryCatch(
      {
        respData = dbGetQuery(db, 
                              paste("SELECT", 
                                    paste0(columns, collapse=','),
                                    " FROM ",
                                    paste0(cte,collapse=" "),
                                    where,
                                    'ORDER BY person_id, item_id;'))
      }, 
      error = function(e)
      {

        if(grepl('(syntax error)|(no such function)', e$message, ignore.case=TRUE))
        {
          # could be that dbplyr cannot translate the r syntax since it contains non sql functions
          # try to read all data (no where clause) and run the select on that
          
          respData <<- dbGetQuery(db, 
                                  paste("SELECT", 
                                        paste0(intersect(union(columns, all.vars(qtpredicate)),
                                                         get_variables(db)$name), collapse=','),
                                        " FROM ",
                                        paste0(cte,collapse=" "),
                                        'ORDER BY person_id, item_id;'))
          
          respData <<- respData[eval_tidy(qtpredicate, data=respData, env=env), columns]
          
        } else
        {
          stop(e)
        }
      },
      finally=NULL)
  }
  return(respData)
}

#### internal ####

# Get respons data including sumscores and design
# 
# @param dataSrc either: data.frame(person_id, item_id, item_score) like what comes out of get_responses,
# or a dexter db handle or an internal `dx_resp_data` object which is trusted for correct booklets and handled efficiently
# @param qtpredicate the quoted predicate, i.e. eval(substitute(quote(predicate)))
# @param extra_columns columns to extract and return on top of the regular ones
# @param extra_design_columns columns to return in the design, can be either in dxItems or dxBooklet_design if dataSrc is a database 
# or they must be part of the input dataframe 
# @param summarised logical, see below
# @param env environment to evaluate the qtpredicate in
#
# @return object of type `dx_resp_data` which is a list containing: 
# `design`: tibble(booklet_id <char or int>, item_id <char or int>, [item position <int> if it can be inferredd])
# `summarised`: logical, copy of parameter
# if summarised == FALSE:
#   `x`: tibble(booklet_id <char or int>, person_id <char>, item_id <char>, item_score <int>, sumScore <int>) 
#   (one row per person-item)
# if summarised == TRUE:
#   `x`: tibble(booklet_id <char or int>, person_id <char>,  sumScore <int>) 
#   (one row per person-booklet)
# @details 
# in case of a qtpredicate that can potentially harm booklets or a data.frame as dataSrc
# this function computes booklets from the data 
get_resp_data = function(dataSrc, qtpredicate=NULL, extra_columns=NULL, extra_design_columns=NULL, summarised=FALSE, env=NULL){
  # gets data and adds sumscores (and booklets if necessary)
  if(is.null(env)) env = caller_env()
  
  ######## resp_data > resp_data ########
  if(inherits(dataSrc,'dx_resp_data'))
  {
    if(!is.null(qtpredicate)) stop("predicates don't work in combination with dx_resp_data object yet")
    
    if(summarised == TRUE & !dataSrc$summarised)
    {
      if(is.null(extra_columns)) 
      {
        dataSrc$x = dataSrc$x %>%
          group_by(.data$booklet_id, .data$person_id) %>% 
          summarise(sumScore = sum(.data$item_score))  %>% 
          ungroup()
      } else
      {
        dataSrc$x = dataSrc$x %>%
          group_by(.data$booklet_id, .data$person_id) %>% 
          mutate(sumScore = sum(.data$item_score)) %>% 
          slice(1) %>%
          ungroup()
      }
      
      dataSrc$summarised = TRUE
    } else if(dataSrc$summarised != summarised)
    {
      stop('cannot unsummarise')
    } 
    
    return(dataSrc)
  }
  
  ######### data.frame or database ##########
  
  bkl_safe = is_bkl_safe(dataSrc, qtpredicate) 

  if(!inherits(dataSrc,'data.frame') || 'booklet_id' %in% colnames(dataSrc))
  {
    # we retrieve the booklet_id column with the data
    columns = c('booklet_id','person_id','item_id','item_score')
  } else
  {
    # no booklet_id column available
    columns = c('person_id','item_id','item_score')
  }
    
  columns = union(columns, extra_columns)

  if(bkl_safe)
  {
    x = get_responses_(dataSrc, qtpredicate = qtpredicate, columns = columns, env = env) 
  } else
  {
    # this also goes for all dataframes since they are never considered safe
    x = get_responses_(dataSrc, qtpredicate = qtpredicate, columns = union(columns, extra_design_columns), env = env) 
  }
  
  
  #### if not bkl_safe, make a booklet_id column
  #### doesn't matter if it is from db or not anymore
  if(!bkl_safe)
  {
    x = x %>%	arrange(.data$person_id, .data$item_id) 
    
    if('booklet_id' %in% colnames(x))
    {
      x$old_booklet_id = factor(x$booklet_id)
    } 
    
    fitem = factor(x$item_id)
    x$iid = fmatch(x$item_id, fitem)
    
    x = x %>% 
      group_by(.data$person_id) %>% 
      mutate(sumScore = sum(.data$item_score), booklet_id = paste0(.data$iid, collapse=' ')) %>% 
      ungroup() 
    
    fbook = factor(x$booklet_id)
    x$booklet_id = fmatch(x$booklet_id, fbook)
    
    if('old_booklet_id' %in% colnames(x))
    {
      # test for surjective mapping, I guess this is enough for calibration but it might not be optimal 
      # since the old booklet can be a subobtimal grouping, that is, there are potentially fewer distinct booklets in the data
      # but if we assume the booklet_id column is at least somewhat serious, it should not be a big problem to err on that side
      if(nlevels(x$old_booklet_id) == nrow(distinct(x, .data$old_booklet_id, .data$booklet_id)))
        x$booklet_id = as.character(x$old_booklet_id)

      # if this also holds then we have a bijective mapping, so it is guaranteed optimal. Should we check for optimal?
      #nrow(distinct(x, .data$booklet_id)) == nrow(distinct(x, .data$booklet_id, .data$old_booklet_id))
    } 


    design = x %>% 
      select(suppressWarnings(one_of(c('booklet_id', 'item_id', 'item_position',extra_design_columns)))) %>%
      distinct(.data$booklet_id, .data$item_id, .keep_all=TRUE)
    
    x = select(x, one_of(c('booklet_id', 'sumScore', columns)))
    
    
    if(summarised == TRUE) {
      x = x %>%
        distinct(.data$booklet_id, .data$person_id, .keep_all=TRUE)
    }
    
  } else 
  {
    ### get design from database ###
    # since it is bkl_safe it must be from a database since data.frames are never trusted
    
    dcol = union(c('booklet_id', 'item_id', 'item_position'), extra_design_columns)
    
    if(length(dcol)==3)
    {
      design = dbGetQuery(dataSrc,'SELECT booklet_id, item_id, item_position FROM dxBooklet_design;')
    } else
    {
      design = dbGetQuery(dataSrc,
                          paste('SELECT',paste0(dcol, collapse=','),
                                'FROM dxBooklet_design INNER JOIN dxItems USING(item_id);'))
    }
    # omit some booklets from design if necessary
    if(!is.null(qtpredicate))
    {
      design = design %>%
        semi_join(x, by='booklet_id')
    }
    
    if(summarised==TRUE)
    {
      if(is.null(extra_columns))
      {
        x = x %>% 
          group_by(.data$person_id, .data$booklet_id) %>% 
          summarise(sumScore = sum(.data$item_score))  %>% 
          ungroup()
      } else {
        x = x %>% 
          select(one_of(c('booklet_id', 'person_id', 'item_score', extra_columns))) %>%
          group_by(.data$person_id, .data$booklet_id) %>% 
          mutate(sumScore = sum(.data$item_score))  %>% 
          slice(1) %>%
          ungroup() %>%
          select(-.data$item_score)
      }
    } else
    {
      x = x %>% 
        group_by(.data$person_id, .data$booklet_id) %>% 
        mutate(sumScore = sum(.data$item_score))  %>% 
        ungroup()
    } 
  }
  
  out = list(x = x, design = arrange(design, .data$booklet_id, .data$item_id), summarised = summarised)
  class(out) = append('dx_resp_data', class(out))
  
  return(out)
}


# filter a dx_resp_data object similar to a dplyr filter statement
# This is less flexible than get_resp_data because
# the predicate here may only contain columns that appear in both the design and the respons data
# on the plus side, if you already have a dx_resp_data object this is much faster than 
# retrieving data from the database
# @parameter respData an object of type dx_resp_data
# @parameter predicate statement (raw/unquoted) to filter data, the predicate may only contain columns that appear in both 
# the design and the data 
# @parameter .recompute_sumscores shall sumScores be recomputed
filter.dx_resp_data = function(respData, predicate, env=NULL, .recompute_sumscores = FALSE)
{
  if(is.null(env)) env = caller_env()
  qtpredicate = eval(substitute(quote(predicate)))
  
  if(.recompute_sumscores & respData$summarised) stop('cannot recompute sumScores on summarised data')
  
  respData$filter = c(respData$filter, as.character(qtpredicate))
  
  # this works but not with the .data pronoun unfortunately
  # respData$design = respData$design[with(respData$design, eval(partial_eval(qtpredicate,env=env))),]
  # respData$x = respData$x[with(respData$x, eval(partial_eval(qtpredicate,env=env))),]
  
  respData$design = respData$design[eval_tidy(qtpredicate, data = respData$design, env = env),]
  respData$x = respData$x[eval_tidy(qtpredicate, data = respData$x, env = env),]
  
  
  if(.recompute_sumscores)
  {
    respData$x = respData$x %>% 
      group_by(.data$person_id, .data$booklet_id) %>% 
      mutate(sumScore = sum(.data$item_score))  %>% 
      ungroup()
  }
  
  return(respData)
}

# filter join for a dx_resp_data object similar to a dplyr semi_join statement
# @parameter respData an object of type dx_resp_data
# @parameter y tibble for semi_join
# @parameter by may only contain columns that are present in the design and in the data(x)
# @parameter .recompute_sumscores shall sumScores be recomputed
semi_join.dx_resp_data = function(respData, y, by=NULL, .recompute_sumscores = FALSE)
{
  if(.recompute_sumscores & respData$summarised) stop('cannot recompute sumScores on summarised data')

  respData$design = semi_join(respData$design, y, by = by)
  respData$x = semi_join(respData$x, y, by = by)
  
  if(.recompute_sumscores)
  {
    respData$x = respData$x %>% 
      group_by(.data$person_id, .data$booklet_id) %>% 
      mutate(sumScore = sum(.data$item_score))  %>% 
      ungroup()
  }
  
  return(respData)
}





