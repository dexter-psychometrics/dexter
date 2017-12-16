

#' Create a database for the 3DC standard setting application
#'
#' Creates an empty database for 3DC standard setting application
#'
#'
#' @param export_name path to a new 3DC database
#' @return a handle to the 3DC sqlite database
#' @details The data driven direct concensus (3DC) method of standard setting  is described in Keuning et. al. (2017).
#' To easily apply this procedure, we advise to use the free digital 3DC application. This application 
#' can be downloaded from the Cito website, see the 
#' \href{http://www.cito.com/our-expertise/implementation/3dc}{3DC application download page}. 
#' The functions create3DC and add_test3DC can be used to produce a standard setting database that 
#' can be imported in the 3DC application.
#' 
#' If you want to apply the 3DC method using paper forms instead, you can use the function plot3DC to generate the forms
#' from the 3DC database.
#' 
#' @references 
#' Keuning J., Straat J.H., Feskens R.C.W. (2017) The Data-Driven Direct Consensus (3DC) Procedure: A New Approach to Standard Setting. 
#' In: Blömeke S., Gustafsson JE. (eds) Standard Setting in Education. 
#' Methodology of Educational Measurement and Assessment. Springer, Cham
#' 
#' @examples
#' \dontrun{
#' library(dplyr)
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' 
#' add_booklet(db, verbAggrData, "aggression")
#' 
#' par = fit_enorm(db) 
#' pv = plausible_values(db, par)
#' mu = mean(pv$PV1)
#' sigma = sd(pv$PV1)
#' 
#' # We'll use the behavior an item depicts as a basis for making the clusters,
#' # thus creating clusters of similar items. 
#' 
#' design = data.frame(item_id = verbAggrProperties$item_id, 
#'    cluster = verbAggrProperties$behavior)
#' 
#' # specify the actual sample for display in the group_leader page
#' 
#' population = get_testscores(db) %>% 
#'                 group_by(test_score) %>% 
#'                 summarise(frequency=n()) 
#' 
#' db3dc = create3DC('test3DC.db')
#' 
#' add_test3DC(db3dc, parms=par, design, mu=mu, sigma=sigma, 
#'             test_id='verbal_aggression', standards='verbally aggressive',
#'             population=population)
#' 
#' #get a preview
#' plot3DC(db3dc)
#' 
#' 
#' dbDisconnect(db3dc)
#' close_project(db)
#'}
#'
create3DC = function(export_name)
{
  if (file.exists(export_name)) file.remove(export_name)
  db3dc = dbConnect(SQLite(), export_name)
  
  dbRunScript(db3dc, '3dc.sql')
  return(db3dc)
}



#' Add a standard setting booklet to a 3DC database
#'
#' See create3DC for more information
#'
#' @param db3dc 3dc database handle
#' @param parms a parameters object produced by fit_enorm
#' @param design a data.frame with columns item_id and cluster and optionally cluster_nbr and item_nbr. See details.
#' @param test_id name/id of the test as it will be shown in the 3DC application.
#' @param standards vector of standards to be set
#' @param mu expected ability in population, used for scaling of the clusters
#' @param sigma expected standard deviation of ability in population, used for scaling of the clusters
#' @param population optionally a data.frame with columns test_score and frequency to use for visualisation in 3DC 
#' application. If NULL, the distribution will be derived from a simulation.
#' @param group_leader Login name of the group leader, default is admin. The default password is always
#' admin, but can be changed in the 3DC application.
#' @param omit the tail probability of the testscores to be omitted. For example, if set to 0.1, 
#' the 10% highest and the 10% lowest scores will be omitted thus restricting the range on which standards can be set.
#' Default is 0.0 (omit nothing)
#'
add_test3DC = function(db3dc, parms, design, test_id, standards, mu, sigma, 
                        population = NULL, group_leader = 'admin', omit = 0.0)
{

  # prevent stupid dplyr warnings
  design$item_id = as.character(design$item_id)
  design$cluster = as.character(design$cluster)
  
  # preserve the order of things if necessary
  if(!'cluster_nbr' %in% names(design)) design$cluster_nbr = match(design$cluster, factor(design$cluster))

  if(!'item_nbr' %in% names(design))
  {
    design = design %>% 
      group_by(.data$cluster_nbr) %>%
      mutate(item_nbr = c(1:n())) %>%
      ungroup()
  }

  m = 300000
  n = nrow(design)
  
  design =  design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  if(parms$input$method=='CML') { b = parms$est$b 
  } else { b = colMeans(parms$est$b) }
  a = parms$inputs$ssIS$item_score
  
  ## Data complete
  x = matrix(0, m, n)
  colnames(x) = design$item_id
  theta = rnorm(m, mu, sigma)

  for (i in 1:n) x[,i] = renorm(b, a, theta, design$first, design$last, i)
  
  cl_scores = as_tibble(sapply(unique(design$cluster), function(cl)
  {
    rowSums(x[,design[design$cluster==cl,]$item_id])
  }, simplify = FALSE, USE.NAMES = TRUE))
  rm(x)
  cl_scores$test_score = rowSums(cl_scores)
  
  max_score = sum((parms$inputs$ssIS %>%  group_by(.data$item_id) %>% 
                    summarise(item_max=max(.data$item_score)) %>%
                    ungroup() %>%
                    inner_join(design,by='item_id'))$item_max)
 
  dbTransaction(db3dc,{
    # add all the test design information to the 3dc database
    dbExecute(db3dc,
              "INSERT INTO Tests (test_id, test_min_score, test_max_score) VALUES(:test_id, 0, CAST(:max_score AS INTEGER));",
              tibble(test_id = test_id, max_score = max_score))
    
    if(is.null(population))
    {
      population = cl_scores %>%
        group_by(.data$test_score) %>%
        summarise(frequnecy = round(n()/30)) # pop of 10.000 seems large enough for histogram
    }  
    population$test_id=test_id
    
    dbExecute(db3dc, 
              "INSERT INTO Population(test_id,test_score,test_score_frequency) 
                VALUES(:test_id, :test_score, :frequency);",
              population)
    
    
    if(!dbExists(db3dc,'SELECT 1 FROM Users WHERE username=?;', group_leader))
      dbExecute(db3dc, 
                "INSERT INTO Users(username, user_password, user_role, user_realname) 
                        VALUES(:uname, :upass, 'group_leader',:uname);",
                  tibble(uname = group_leader, upass='$5$rounds=110000$XSVXY0auWBdJRaqa$Mp4aIPQoP8281M9FHYfwCWfvxpkht8TVVwxa2XWAhq7')
                )
    
    dbExecute(db3dc, "INSERT INTO group_leader_test_assignment(username, test_id) VALUES(:uname, :test);",
              tibble(uname = group_leader, test = test_id))
    
    
    # reverse the cluster nbr to order them correctly on the screen in 3DC program (historic reasons)
    design$cluster_nbr =  1 + max(design$cluster_nbr) - design$cluster_nbr
    design$test_id = test_id
    
    dbExecute(db3dc, 
              "INSERT INTO Clusters (test_id, cluster_nbr, cluster_name) 
                VALUES(:test_id, CAST(:cluster_nbr AS INTEGER), :cluster);",
              design[,c('test_id', 'cluster_nbr', 'cluster')] %>% 
                group_by(.data$test_id, .data$cluster_nbr, .data$cluster) %>%
                slice(1))
  
    dbExecute(db3dc,
              "INSERT INTO Standards(test_id, standard_nbr, standard_name) 
                VALUES(:test_id, CAST(:nbr AS INTEGER), :name);",
              tibble(test_id = test_id, nbr = seq_along(standards), name = standards))
    
    included_scores =
      cl_scores  %>% 
      arrange(.data$test_score) %>% 
      group_by(.data$test_score) %>%
      summarise(N=n()) %>% 
      filter((cumsum(.data$N) >= omit*m) & ((m - cumsum(.data$N)) >= omit*m))  %>%
      select(.data$test_score)
    if(nrow(included_scores) == 0) stop('all scores are exluded')
    
    # estimates C+|X+
    cl_scores  = cl_scores %>% 
      inner_join(included_scores, by='test_score') %>%
      group_by(.data$test_score) %>%
      do(as_tibble(t(floor(apply(.[,unique(design$cluster)], 2, mean))))) %>%
      ungroup()
  
    apply(unique(design[,c('cluster','cluster_nbr')]), 1, function(cl)
    {
      by(cl_scores, cl_scores[,cl['cluster']], function(score)
      {
        dbExecute(db3dc,
                  "INSERT INTO Cluster_scores(test_id, cluster_nbr, cluster_score,  test_score_est)
                  VALUES(:test_id, :cluster_nbr, :cl_score, :test_score_est);",
                  tibble(test_id = test_id, cluster_nbr = cl['cluster_nbr'],
                             cl_score = unique(score[[cl['cluster']]]),
                             test_score_est = min(score['test_score'])))
      })
    })
    
    for(itm in unique(design$item_id))
    {
      if(!dbExists(db3dc, 'SELECT 1 FROM Items WHERE item_id=?;', itm))
      {
        dbExecute(db3dc, 'INSERT INTO Items(item_id, item_name) VALUES(:itm, :itm);', tibble(itm=itm))
      }
    }
    
    dbExecute(db3dc,
              "INSERT INTO Cluster_design (test_id, cluster_nbr, item_nbr, item_id) 
                VALUES(:test_id, :cluster_nbr, :item_nbr, :item_id);",
                design[,c('test_id', 'cluster_nbr', 'item_nbr', 'item_id')])
    
  })
  invisible(NULL)
}


#' Show 3DC plots
#'
#' Show 3DC plots as used in 3DC standard setting application
#'
#'
#' @param db3dc 3dc database handle
#' @param test_id optionally, a vector of test_id's. If omitted, plots for all tests will be shown
#' @param ... further arguments to plot. Some of them have useful defaults
#' 
plot3DC = function(db3dc, test_id=NULL, ...)
{
  where = ''
  if(!is.null(test_id)) where = paste0('WHERE test_id IN(',paste0(sql_quote(test_id), collapse=','),')')
  
  plot.args = list(...)
  
  dm = dbGetQuery(db3dc,paste('SELECT test_id, test_max_score, COUNT(*) AS N 
                                FROM Tests INNER JOIN Clusters USING(test_id)',
                                where,
                              'GROUP BY test_id, test_max_score;'))
             
  
  plt = dbGetQuery(db3dc, paste('SELECT test_id, cluster_name, cluster_nbr, cluster_score,  test_score_est 
                                  FROM Cluster_scores 
                                    INNER JOIN Clusters USING(test_id, cluster_nbr)',
                                  where,
                                  'ORDER BY test_id, cluster_nbr DESC, cluster_score;'))
  
  clusters = dbGetQuery(db3dc, 'SELECT test_id, cluster_name, cluster_nbr FROM Clusters;')

  
  plt %>% 
    group_by(.data$test_id) %>%
    do({
      test = .$test_id[1]
      maxSc = dm[dm$test_id == test,]$test_max_score 
      N = dm[dm$test_id == test,]$N
      default.args = list(
        xaxp = c(0,maxSc,maxSc),
        xlim = c(0,maxSc),
        ylim = c(0,N+.3),
        xlab = 'Score',
        ylab = 'Cluster',
        main = test,
        axes = FALSE
      )
      do.call(plot,
              merge_arglists(plot.args,
                             override = list(x = c(NA,NA), y = c(maxSc,N)),
                             default = default.args))
      axis(1,at=c(0:maxSc))
      #axis(2,at=c(1:N))

      for(i in c(1:N))
      {
        graphics::lines(c(0,maxSc),c(i,i))
        graphics::text(1,i,clusters[clusters$test_id==test & clusters$cluster_nbr==i,]$cluster_name,adj=c(0,-1))
      }
      apply(.,1, function(x)
      {
        graphics::text(as.integer(x['test_score_est']), as.integer(x['cluster_nbr']),x['cluster_score'], pos=1)
        graphics::points(as.integer(x['test_score_est']), as.integer(x['cluster_nbr']),pch=20)
      })
      data.frame()
    })
  invisible(NULL)
}

