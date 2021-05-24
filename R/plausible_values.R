

#' Draw plausible values
#'
#' Draws plausible values based on test scores
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by function \code{fit_enorm} containing parameter estimates. If parms are provided, item parameters are considered known. 
#' If parms = NULL, plausible values are marginalized over the posterior distribution of the item parameters and uncertainty of the item parameters is taken into account.
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param covariates name or a vector of names of the variables to group the populations used to improve the prior.
#' A covariate must be a discrete person property (e.g. not a float) that indicates nominal categories, e.g. gender or school.
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param nPV Number of plausible values to draw per person.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler, this specifies the use of a particular sample of item parameters used to generate the plausible value(s). If NULL, the posterior means are used. If outside range, the last iteration will be used.
#' @param prior.dist use a normal prior or a mixture of two normals
#' recognised automatically), 
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible values are generated per person (TRUE) or per booklet (FALSE)
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
                            merge_within_persons=FALSE)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  prior.dist = match.arg(prior.dist)
  check_dataSrc(dataSrc)
  check_num(nPV, .length=1, .min=1)
  
  df_format(
    mutate_if(
      plausible_values_(dataSrc, parms, qtpredicate=qtpredicate, covariates=covariates, nPV=nPV, 
                      use_draw=use_draw, env=env,prior.dist = prior.dist ,
                      merge_within_persons=merge_within_persons), 
    is.factor, as.character))
}

# to~do: ignore covariate when (some) groups contain to few, <5 say, persons. Add warning.
# what if these are 4 persons with score 0 on an easy test?
# would, in general, the proper way to deal with the pathological case be to add a dummy covariate
# based on characteristics of scoretab? (per booklet and per user covariate of course) 

plausible_values_ = function(dataSrc, parms=NULL, qtpredicate=NULL, covariates=NULL, nPV=1, use_draw=NULL, 
                             env=NULL, prior.dist = c("normal", "mixture"),
                             merge_within_persons=merge_within_persons)
{
  if(is.null(env)) env = caller_env()
  from = Gibbs.settings$from.pv
  step = Gibbs.settings$step.pv # burnin and thinning for pvs
  nIter.enorm = from + step*(nPV-1) # nr. of posterior samples of item parameters needed
  
  prior.dist = match.arg(prior.dist)
  
  pb = get_prog_bar(nsteps=if(is.null(parms)) 100 else 120, 
                    retrieve_data = is_db(dataSrc))
  on.exit({close_prog_bar()})

  if(is.null(parms))
  {
    respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env)
    pb$open_sub_bar(20)
    parms = fit_enorm_(respData, method = 'Bayes', nDraws = nIter.enorm) 
    respData = get_resp_data(respData, summarised=TRUE, extra_columns=covariates, 
                             protect_x=!is_db(dataSrc))
    pb$open_sub_bar(100)

  } else
  {
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
  
  
  # to do: use simplify_parms
  if(inherits(parms,'data.frame'))
  {
    fl = parms %>%
      mutate(rn=row_number()) %>%
      group_by(.data$item_id) %>%
      summarize(first=as.integer(min(.data$rn)), last=as.integer(max(.data$rn))) %>%
      ungroup()
    a = parms$item_score
    b = parms$b
  } else
  {
    fl = parms$inputs$ssI[c('item_id','first','last')]
    a = parms$inputs$ssIS$item_score
    b = parms$est$b
    if(parms$inputs$method == 'Bayes' )
    {
      if(!is.null(use_draw))
      {
        b = b[use_draw,]
      } else if (nrow(b)<nIter.enorm ) 
      {
        stop(paste("To produce", nPV, "plausible values, use at least", nIter.enorm, "iterations in fit_enorm" ))
      } 
    }  
  }
  x = respData$x
  
  # join design with the params
  # these can have different levels
  design = suppressWarnings(respData$design %>%
    inner_join(fl, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first))

  if(!is.null(covariates))
  {
    group_number = (function(){i = 0L; function() i <<- i+1L })()
    x = x %>% 
      group_by_at(covariates) %>%
      mutate(pop__ = group_number()) %>%
      ungroup() 
  } else
  {
    # niet varierende pop toevoegen maakt code in pv eenvoudiger
    x$pop__ = 1L
  }
  design = split(design, design$booklet_id, drop=TRUE)
  
  y = pv(select(x, .data$booklet_id, .data$person_id, .data$booklet_score, pop = .data$pop__),
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

