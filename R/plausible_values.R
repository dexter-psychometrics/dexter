
#' Draw plausible values
#'
#' Draws plausible values based on test scores
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by function \code{fit_enorm} containing parameter estimates. If parms are provided, item parameters are considered known. 
#' If parms = NULL, they will be estimated Bayesianly.
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param covariates name or a vector of names of the variables to group the populations used to improve the prior.
#' A covariate must be a discrete person property (e.g. not a float) that indicates nominal categories, e.g. gender or school.
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param nPV Number of plausible values to draw per person.
#' @param parms_draw when the item parameters are estimated Bayesianly (see: \code{\link{fit_enorm}}), 
#' parms_draw specifies whether to use a sample(a different item parameter draw for each plausible values draw) or the posterior mean
#' of the item draws. Alternatively, it can be an integer specifying a specific draw. Ignored when parms is not estimated Bayesianly.
#' @param prior_dist use a normal prior or a mixture of two normals
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible values are generated per person (TRUE) or per booklet (FALSE)
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPV plausible values
#' named PV1...PVn.
#' 
#' @details
#' 
#' When The item parameters are estimated using \code{fit_enorm(..., method='Bayes')} and parms_draw = 'sample', the uncertainty 
#' of the item parameters estimates is taken into account when drawing multiple plausible values. 
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
plausible_values_old = function(dataSrc, parms=NULL, predicate=NULL, covariates=NULL, 
                            nPV=1, 
                            parms_draw = c('sample','average'), 
                            prior_dist = c("normal", "mixture"),
                            merge_within_persons = FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  prior_dist = match.arg(prior_dist)
  check_dataSrc(dataSrc)
  check_num(nPV, .length=1, .min=1)
  
  plausible_values_old_(dataSrc, parms, qtpredicate=qtpredicate, covariates=covariates, nPV=nPV, 
                    parms_draw = parms_draw, env=env,prior_dist = prior_dist ,
                    merge_within_persons=merge_within_persons) %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}

# to~do: ignore covariate when (some) groups contain to few, <5 say, persons. Add warning.
# what if these are 4 persons with score 0 on an easy test?
# would, in general, the proper way to deal with the pathological case be to add a dummy covariate
# based on characteristics of scoretab? (per booklet and per user covariate of course) 

plausible_values_old_ = function(dataSrc, parms=NULL, qtpredicate=NULL, covariates=NULL, nPV=1, parms_draw = c('sample','average'), 
                             env=NULL, prior_dist = c("normal", "mixture"),
                             merge_within_persons=merge_within_persons)
{
  if(is.null(env)) env = caller_env()
  
  if(is.numeric(parms_draw)) parms_draw = as.integer(parms_draw)
  else parms_draw = match.arg(parms_draw)
  
  prior_dist = match.arg(prior_dist)
  
  pv_from = Gibbs.settings$from.pv 
  pv_step = Gibbs.settings$step.pv
  niter_req = pv_from + pv_step*(nPV-1) 
  
  pb = get_prog_bar(nsteps=if(is.null(parms)) 120 else 100, 
                    retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  if(is.null(parms))
  {
    respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env)
    pb$new_area(20)
    parms = fit_enorm_(respData, method = 'Bayes', nDraws = niter_req) 
    
    respData = get_resp_data(respData, summarised=TRUE, extra_columns=covariates, 
                             protect_x=!is_db(dataSrc))
    pb$new_area(100)
    
  } else
  {
    # to do: can simplify parms be done here?
    if(inherits(parms,'data.frame'))
    {
      parms = transform.df.parms(parms,'b', TRUE)
      pcheck = parms[,c('item_id','item_score')]
    } else
    {
      pcheck = parms$inputs$ssIS[,c('item_id','item_score')]
    }
    
    respData = get_resp_data(dataSrc, qtpredicate, summarised=TRUE, 
                             extra_columns=covariates,env=env, 
                             parms_check=pcheck,
                             merge_within_persons=FALSE)
  }
  
  parms = simplify_parms(parms, draw=parms_draw)
  
  if(parms_draw == 'sample' && parms$method != 'CML')
  {
    if(nrow(parms$b) < niter_req )
      stop_(paste("To produce", nPV, "plausible values, use at least", niter_req, "iterations in fit_enorm" ))
  }
  if(!is.null(covariates))
  {
    group_number = (function(){i = 0L; function() i <<- i+1L })()
    respData$x = respData$x %>% 
      group_by_at(covariates) %>%
      mutate(pop__ = group_number()) %>%
      ungroup() 
  } else
  {
    # niet varierende pop toevoegen maakt code in pv eenvoudiger
    respData$x$pop__ = 1L
  }
  
  # join design with the params
  # these can have different levels
  design = suppressWarnings(respData$design %>%
                              inner_join(parms$items, by='item_id') %>% 
                              arrange(.data$booklet_id, .data$first))
  
  design = split(design, design$booklet_id, drop=TRUE)
  
  y = pv(select(respData$x, 'booklet_id', 'person_id', 'booklet_score', pop = 'pop__'),
         design, parms$b, parms$a, nPV, from = pv_from, by = pv_step, prior.dist = prior_dist)
  
  colnames(y) = c('booklet_id','person_id','booklet_score',paste0('PV',1:nPV))
  
  if(is.null(covariates))
  {
    y
  } else
  {
    # added unique so that booklet_id can be used as a covariate
    inner_join(respData$x[,unique(c('booklet_id','person_id',covariates))], y, by=c('booklet_id','person_id') )
  }
}

