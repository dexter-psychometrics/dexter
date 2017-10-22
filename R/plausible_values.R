

##########################################
#' Draw plausible values
#'
#' Draws plausible values based on sum scores
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates. If parms is given the function provides plausible values conditional on the 
#' item paramaters. These are considered known. If parameter estimates are not given, the user is given 
#' plausible value marginalized over the posterior distribution of the item parameters. 
#' In plain words, this means that the uncertainty of the item parameters is taken into account.
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param covariates name or a vector of names of the variables to group the population, used to update the prior.
#' A covariate must be a discrete person covariate (e.g. not a float) that indicates nominal categories, e.g. gender or school
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param nPV Number of plausible values to draw per person.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler (this is 
#' recognised automatically), the number of the random draw (iteration) to use 
#' in generating the PV. If NULL, all draws will be averaged; that is, the posterior means are used for the item parameters.
#' If outside range, the last iteration will be used.   
#' @return A data.frame with columns booklet_id, person_id, sumScore and nPV plausible values
#' named PV1...PVn.
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
#' plot(ecdf(pv_M$PV1), 
#'    main="Do: males versus females", xlab="Ability", col="red")
#' lines(ecdf(pv_F$PV1), col="green")
#' legend(-2.2,0.9, c("female", "male") , 
#'    lty=1, col=c('green', 'red'), bty='n', cex=.75)
#'
#' pv_M=plausible_values(db,f,(mode=="Want")&(gender=="Male"))
#' pv_F=plausible_values(db,f,(mode=="Want")&(gender=="Female"))
#'
#' plot(ecdf(pv_M$PV1), 
#'    main="Want: males versus females", xlab=" Ability", col="red")
#' lines(ecdf(pv_F$PV1),col="green")
#' legend(-2.2,0.9, c("female", "male") , 
#'    lty=1, col=c('green', 'red'), bty='n', cex=.75)
#'    
#' close_project(db)    
#' }
#' 
plausible_values = function(dataSrc, parms=NULL, predicate=NULL, covariates=NULL, nPV=1, use_draw=NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  plausible_values_(dataSrc, parms, qtpredicate=qtpredicate, covariates=covariates, nPV=nPV, use_draw=use_draw, env=caller_env())
}

plausible_values_ = function(dataSrc, parms, qtpredicate=NULL, covariates=NULL, nPV=1, use_draw=NULL, env=NULL)
{
  if(is.null(env)) env = caller_env()

  parms_given = !is.null(parms)

  if (parms_given)
  {
    r = get_resp_data(dataSrc, qtpredicate, summarised=TRUE, extra_columns=covariates,env=env)
  } else
  {
    r = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env)
    parms = fit_enorm_(r, method='Bayes', nIterations=(4+nPV*5))
    r = get_resp_data(r, summarised=TRUE, extra_columns=covariates)
  } 
 
  x = r$x
  if(nrow(x) == 0) stop('no data to analyse')
  
  # join design with the params
  design = r$design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  if (parms_given)
  {
    if(parms$input$method=='CML') 
    {
      b = parms$est$b
      a = parms$inputs$ssIS$item_score
    } else
    { 
      a = parms$est$a
      if(is.null(use_draw)) 
      {
        b = colMeans(parms$est$b)  
      } else 
      {
        if (use_draw %in% 1:nrow(parms$est$b)) 
        {
          b = parms$est$b[use_draw,]  
        } else 
        {
          b = parms$est$b[nrow(parms$est$b),]
        }
      }   
    }
  }else
  {
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
  }
    
    
  if(!is.null(covariates))
  {
    group_number = (function(){i = 0; function() i <<- i+1 })()
    
    x = x %>% 
      group_by_at(covariates) %>%
      mutate(pop = group_number()) %>%
      ungroup()
  } else
  {
    # niet varierende pop toevoegen maakt code in pv eenvoudiger
    x$pop = 1
  }
  # design as list makes pv faster
  bkl = unique(design$booklet_id)
  design = lapply(bkl, function(x){
    design %>% filter(.data$booklet_id==x) %>% select(.data$first, .data$last) %>% arrange(.data$first)
  })
  names(design) = as.character(bkl)
  
  # arranging probably makes pv faster but it is not required
  y = x %>% arrange(.data$booklet_id, .data$pop) %>% pv(design, b, a, nPV)
  
  colnames(y) = c('booklet_id','person_id','sumScore',paste0('PV',1:nPV))
  if(is.null(covariates))
  {
    return(y)
  } else
  {
    return(inner_join(y, x[,c('booklet_id','person_id',covariates)], by=c('booklet_id','person_id')))
  }
}

