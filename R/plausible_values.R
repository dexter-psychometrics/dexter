

##########################################
#' Draw plausible values
#'
#' Draws plausible values based on test scores
#'
#'
#' @param dataSrc Data source: a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates. If parms is given the function provides plausible values conditional on the 
#' item parameters; i.e., these are considered known and might be based on a different data set. 
#' If parms = NULL, the user is given plausible values marginalized over the posterior distribution of the item parameters. 
#' In plain words, this means that the uncertainty of the item parameters is taken into account.
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param covariates name or a vector of names of the variables to group the populations used to improve the prior.
#' A covariate must be a discrete person property (e.g. not a float) that indicates nominal categories, e.g. gender or school
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param nPV Number of plausible values to draw per person.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler (this is 
#' @param prior.dist use a normal prior or a mixture of two normals
#' recognised automatically), the number of the random draw (iteration) to use 
#' in generating the PV. If NULL, all draws will be averaged; that is, the posterior means are used for the item parameters. If outside range, the last iteration will be used.
#' @param merge_within_person for persons who were administered multiple booklets, 
#' whether to provide just one plausible value (TRUE) or one per booklet(FALSE)
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPV plausible values
#' named PV1...PVn.
#' 
#' @references 
#' Marsman, M., Maris, G., Bechger, T. M., and Glas, C.A.C. (2016) What can we learn from plausible values? 
#' Psychometrika. 2016; 81: 274-289. See also the vignette.
#' 
#' @examples
#' db = start_new_project(verbAggrRules, ":memory:", 
#'    person_properties=list(gender="<unknown>"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' 
#' f=fit_enorm(db)
#' pv_M=plausible_values(db,f,(mode=="Do")&(gender=="Male"))
#' pv_F=plausible_values(db,f,(mode=="Do")&(gender=="Female"))
#' 
#' par(mfrow=c(1,2))
#' 
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
#' 
plausible_values = function(dataSrc, parms=NULL, predicate=NULL, covariates=NULL, 
                            nPV=1, use_draw=NULL, prior.dist = c("normal", "mixture"),
                            merge_within_person=FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  prior.dist = match.arg(prior.dist)
  check_dataSrc(dataSrc)
  check_num(nPV, .length=1, .min=1)
  
  plausible_values_(dataSrc, parms, qtpredicate=qtpredicate, covariates=covariates, nPV=nPV, 
                    use_draw=use_draw, env=env,prior.dist = prior.dist ,
                    merge_within_person=merge_within_person) %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
}
# to do: check factor warnings when parms has different items than dataSrc
# to do: ignore covariate when (some) groups contain to few, <5 say, persons. Add warning.
# what if these are 4 persons with score 0 on an easy test?
# would, in general, the proper way to deal with the pathological case be to add a dummy covariate
# based on characteristics of scoretab? (per booklet and per user covariate of course) 

# use_b_matrix makes the function act as if parms=null when parms is actually bayesian
plausible_values_ = function(dataSrc, parms=NULL, qtpredicate=NULL, covariates=NULL, nPV=1, use_draw=NULL, 
                             env=NULL, use_b_matrix=FALSE, prior.dist = c("normal", "mixture"),
                             merge_within_person=merge_within_person)
{
  if(is.null(env)) env = caller_env()
  from = 20 ; step = 5
  nIter_enorm = from + step*(nPV-1)
  prior.dist = match.arg(prior.dist)

  parms_given = !is.null(parms)
  if (parms_given)
  {
    respData = get_resp_data(dataSrc, qtpredicate, summarised=TRUE, 
                             extra_columns=covariates,env=env, 
                             parms_check=parms$inputs$ssIS[,c('item_id','item_score')],
                             merge_within_person=FALSE)
  } else
  {
    respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env)
    parms = fit_enorm_(respData, method = 'Bayes', nIterations = nIter_enorm) 
    respData = get_resp_data(respData, summarised=TRUE, extra_columns=covariates, protect_x=!inherits(dataSrc,'DBIConnection'))
  } 

  x = respData$x
  if(nrow(x) == 0) stop('No data to analyse')
  
  # join design with the params
  # these can have different levels
  design = suppressWarnings(respData$design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first))
  
  if(any(is.na(design$first))) stop('Some of your items are without parameters')
  
  if (parms_given && !use_b_matrix)
  {
    if(parms$inputs$method=='CML') 
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
        b = parms$est$b[min(use_draw,nrow(parms$est$b)),]  
      }   
    }
  }else
  {
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
    
    if (nrow(b)<nIter_enorm) stop(paste("To produce", as.character(nPV), "plausible values. Do at least", 
                                      as.character(nIter_enorm), "iterations of fit_enorm" ))
  }
  
  
  if(!is.null(covariates))
  {
    group_number = (function(){i = 0L; function() i <<- i+1L })()
    
    x = x %>% 
      group_by_at(covariates) %>%
      mutate(pop = group_number()) %>%
      ungroup() %>%
      arrange(.data$booklet_id, .data$pop) 
    # to do: arrange by sumscore in anon (or if(unsorted...)) so it is explicit that sorted pv can be used
    # also possible to make a scoretab at the start
    # to do: check how long this sort takes, might be more efficient in c
    
  } else
  {
    if(is.unsorted(x$booklet_id))
      x = arrange(x, .data$booklet_id)
    
    # niet varierende pop toevoegen maakt code in pv eenvoudiger
    x$pop = 1L
  }
  # design as list makes pv faster
  design = split(design, design$booklet_id)
  
  # from and by are burnin and thinning
  y = pv(x[,c('booklet_id','person_id','booklet_score','pop')], 
         design, b, a, nPV, from = from, by = step, prior.dist=prior.dist)
  
  colnames(y) = c('booklet_id','person_id','booklet_score',paste0('PV',1:nPV))
  
  if(is.null(covariates))
  {
    y
  } else
  {
    # added unique so that booklet_id can be used as a covariate
    inner_join(x[,unique(c('booklet_id','person_id',covariates))], y, by=c('booklet_id','person_id') )
  }
}

