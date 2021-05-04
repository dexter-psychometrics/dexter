

#' Standard setting
#'
#' Set performance standards on one or more test forms using the data driven direct consensus (3DC) method
#'
#' @param parms parameters object returned from fit_enorm
#' @param design a data.frame with columns `cluster_id`, `item_id` and optionally `booklet_id`
#' @param x an object containing parameters for the 3DC standard setting procedure
#' @param object an object containing parameters for the 3DC standard setting procedure
#' @param booklet_id which test form to plot
#' @param ... ignored
#' Optionally you can include a column `booklet_id` to specify multiple test forms for standard setting and/or
#' columns `cluster_nbr` and `item_nbr` to specify ordering of clusters and items in the forms and application.
#' @return 
#' an object of type `sts_par`
#'
#' @details The data driven direct consensus (3DC) method of standard setting was invented by Gunter Maris and described in Keuning et. al. (2017).
#' To easily apply this procedure, we advise to use the free digital 3DC application. This application 
#' can be downloaded from the Cito website, see the 
#' \href{https://www.cito.com/our-expertise/implementation/3dc}{3DC application download page}. 
#' If you want to apply the 3DC method using paper forms instead, you can use the function plot3DC to generate the forms
#' from the 3DC database.
#' 
#' Although the 3DC method is used as explained in Keuning et. al., the method we use for computing the forms is a simple
#' maximum likelihood scaling from an IRT model, described in Moe and Verhelst (2017)
#' 
#' @references 
#' Keuning J., Straat J.H., Feskens R.C.W. (2017) The Data-Driven Direct Consensus (3DC) Procedure: A New Approach to Standard Setting. 
#' In: Blomeke S., Gustafsson JE. (eds) Standard Setting in Education. 
#' Methodology of Educational Measurement and Assessment. Springer, Cham
#' 
#' Moe E., Verhelst N. (2017) Setting Standards for Multistage Tests of Norwegian for Adult Immigrants 
#' In: Blomeke S., Gustafsson JE. (eds) Standard Setting in Education. 
#' Methodology of Educational Measurement and Assessment. Springer, Cham
#' 
#' @seealso how to make a database for the 3DC standard setting application: \code{\link{standards_db}}
#' 
#' @examples
#' 
#' library(dplyr)
#' db = start_new_project(verbAggrRules, ":memory:")
#'             
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' 
#' design = get_items(db) %>%
#'   rename(cluster_id='behavior')
#' 
#' f = fit_enorm(db)
#' 
#' sts_par = standards_3dc(f, design)
#' 
#' plot(sts_par)
#' 
#' 
#' # db_sts = standards_db(sts_par,'test.db',c('mildly aggressive','dangerously aggressive'))
standards_3dc = function(parms, design)
{
  
  check_df(design, c('cluster_id','item_id'))
  design = mutate_if(design, is.factor, as.character)
  
  if(!'booklet_id' %in% colnames(design))
    design$booklet_id = '3DC'

  if(!'cluster_nbr' %in% colnames(design)) 
    design = design %>% 
      group_by(.data$booklet_id) %>%
      mutate(cluster_nbr = dense_rank(.data$cluster_id))

  if(!'item_nbr' %in% names(design))
    design = design %>% 
      group_by(.data$booklet_id, .data$cluster_nbr) %>%
      mutate(item_nbr = dense_rank(.data$item_id)) %>%
      ungroup()


  if(length(setdiff(design$item_id, as.character(parms$inputs$ssI$item_id)))>0)
  {
    warning('The following items are present in design but do not have parameters')
    print(setdiff(design$item_id, as.character(parms$inputs$ssI$item_id)))
    stop('unknown items in design')
  }
  est = lapply(split(design, design$booklet_id), function(tds)
  {
      es = expected_score(parms, items=tds$item_id)
      
      select(tds, booklet_id='cluster_id', .data$item_id) %>%
        ability_tables(parms, design = ., standard_errors=FALSE) %>%
        rename(cluster_id='booklet_id', cluster_score='booklet_score') %>%
        mutate(booklet_score = es(.data$theta)) %>%
        inner_join(distinct(tds, .data$cluster_nbr, .data$cluster_id), by='cluster_id')
  })
  
  
  out = list(design=design, est=est)
  class(out) = append('sts_par', class(out))
  out
}

#' Export a standard setting database for use by the free 3DC application
#' 
#' This function creates an export (an sqlite database file) which can be used by the 3DC application. This is a free application with which
#' a standard setting session can be facilitated through a LAN network using the Chrome browser.
#' The 3DC application can be downloaded from \href{https://www.cito.com/our-expertise/implementation/3dc}{3DC application download page}
#' 
#' 
#' @param par.sts an object containing parameters for the 3DC standard setting procedure produced by
#' \code{\link{standards_3dc}}
#' @param file_name name of the exported database file
#' @param standards vector of 1 or more standards. In case there are multiple test forms and
#' they should use different performance standards, a list of such vectors. 
#' The names of this list should correspond to the names of the testforms
#' @param population optional, a data.frame with three columns: `booklet_id`,`booklet_score`,`n` (where n is a count)
#' @param group_leader login name of the group leader. The login password will always be `admin` 
#' but can be changed in the 3DC application
#' 
standards_db = function(par.sts, file_name, standards, population=NULL, group_leader = 'admin')
{
  if (file.exists(file_name)) 
    if(!file.remove(file_name))
      stop('file already exists and cannot be removed')
  
  booklets = names(par.sts$est)
  if(!is.list(standards))
  {
    standards = rep(list(standards), length(booklets))
    names(standards) = booklets
  }
  if(length(standards) != length(booklets))
      stop('expected standard to be a list of length ', length(booklets), ' found ', length(standards))
  

  db3dc = dbConnect(SQLite(), file_name)
  
  dbRunScript(db3dc, '3dc.sql')
  
  dbTransaction(db3dc,
  {
    dbExecute(db3dc, 
              "INSERT INTO Users(username, user_password, user_role, user_realname) 
                    VALUES(:uname, :upass, 'group_leader',:uname);",
                tibble(uname = group_leader, upass='$5$rounds=110000$XSVXY0auWBdJRaqa$Mp4aIPQoP8281M9FHYfwCWfvxpkht8TVVwxa2XWAhq7'))
    
    dbExecute(db3dc,
                "INSERT INTO Tests (test_id, test_min_score, test_max_score) VALUES(:test_id, 0, CAST(:max_score AS INTEGER));",
                tibble(test_id = names(par.sts$est), max_score = sapply(par.sts$est, function(x) max(x$booklet_score))))
    
    dbExecute(db3dc, "INSERT INTO group_leader_test_assignment(username, test_id) VALUES(:uname, :booklet_id);",
                tibble(uname = group_leader, booklet_id=booklets))
    
    for(booklet in names(standards))
    {
      dbExecute(db3dc,
              "INSERT INTO Standards(test_id, standard_nbr, standard_name) 
                VALUES(:booklet_id, CAST(:nbr AS INTEGER), :standard_name);",
              tibble(booklet_id=booklet, 
                     nbr = seq_along(standards[[booklet]]), standard_name = standards[[booklet]]))
    }
    
    dbExecute(db3dc, 'INSERT INTO Items(item_id, item_name) VALUES(:item_id, :item_id);', 
              distinct(par.sts$design, .data$item_id))
    
    lapply(split(par.sts$design, par.sts$design$booklet_id), function(tds)
    {
      dbExecute(db3dc, 
                "INSERT INTO Clusters (test_id, cluster_nbr, cluster_name) 
                VALUES(:booklet_id, CAST(:cluster_nbr AS INTEGER), :cluster_id);",
                distinct(tds, .data$booklet_id, .data$cluster_nbr, .data$cluster_id))
      
      dbExecute(db3dc,
                "INSERT INTO Cluster_design (test_id, cluster_nbr, item_nbr, item_id) 
                VALUES(:booklet_id, :cluster_nbr, :item_nbr, :item_id);",
                tds[,c('booklet_id', 'cluster_nbr', 'item_nbr', 'item_id')])
      
      x = par.sts$est[[tds$booklet_id[1]]]
      dbExecute(db3dc,
                "INSERT INTO Cluster_scores(test_id, cluster_nbr, cluster_score,  test_score_est)
                  VALUES(:booklet_id, :cluster_nbr, :cluster_score, :booklet_score);",
                tibble(booklet_id=tds$booklet_id[1],
                       cluster_nbr=x$cluster_nbr, cluster_score=x$cluster_score, booklet_score=x$booklet_score))
    })
    

    if(!is.null(population))
    {
      check_df(population, c('booklet_id', 'booklet_score', 'n'))
      #to~do: possible pop from prms?
      if(length(intersect(booklets,population$booklet_id)) < length(booklets))
        stop("one or more booklet_id's in your population do not exist in your sts parameters")
      dbExecute(db3dc, 
              "INSERT INTO Population(test_id,test_score,test_score_frequency) 
                VALUES(:booklet_id, :booklet_score, :n);",
                select(population, .data$booklet_id, .data$booklet_score, .data$n))
    } 
  }, on_error=function(e){
    dbDisconnect(db3dc); file.remove(file_name); stop(e)})
    
  return(db3dc)
}

#'@rdname standards_3dc
coef.sts_par = function(object, ...)
{
  bind_rows(object$est, .id='booklet_id') %>%
    inner_join(distinct(object$design, .data$booklet_id, .data$cluster_nbr, .data$cluster_id), 
               by=c('booklet_id','cluster_nbr','cluster_id')) %>%
    arrange(.data$booklet_id, .data$cluster_nbr, .data$cluster_score) %>%
    df_format()
}

#'@rdname standards_3dc
plot.sts_par = function(x, booklet_id=NULL, ...)
{
  plot.args = list(...)
  if(is.null(booklet_id))
     booklet_id = names(x$est)
     
  for(tst in booklet_id)
  {

    max_score = max(x$est[[tst]]$booklet_score)
    n_cluster = max(x$est[[tst]]$cluster_nbr)
    
    default.args = list(
      xaxp = c(0,max_score,max_score),
      xlim = c(0,max_score),
      ylim = c(0,n_cluster+.3),
      xlab = 'Test score',
      ylab = 'Cluster',
      main = tst,
      axes = FALSE
    )
    
    do.call(plot,
            merge_arglists(plot.args,
                           override = list(x = c(0,max_score), y = c(0, n_cluster), type='n'),
                           default = default.args))
    
    cnames = filter(x$design, .data$booklet_id==tst) %>%
      distinct(.data$cluster_nbr, .data$cluster_id) %>%
      arrange(.data$cluster_nbr) %>%
      pull('cluster_id')
    
    axis(1,at=0:max_score)

    for(i in c(1:n_cluster))
    {
      lines(c(0,max_score),c(i,i))
      text(1,i+.1,cnames[i],adj=c(0,0))
      p = filter(x$est[[tst]], .data$cluster_nbr==i)
      text(p$booklet_score, i, p$cluster_score, pos=1)
      points(p$booklet_score, rep(i,nrow(p)))
    }
  }
}

