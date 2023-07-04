

pv_design = function(design, a)
{
  design = arrange(design, .data$booklet_id, .data$first) |>
    mutate(first_c = as.integer(.data$first-1L),
           last_c = as.integer(.data$last-1L))
  
  booklets = design %>%
    group_by(.data$booklet_id) %>%
    summarise(nit=n(), max_a = max(a[.data$last])) %>%
    ungroup() %>%
    arrange(.data$booklet_id)
  
  list(first_c = design$first_c, last_c = design$last_c, 
       bk_cnit = c(0L,cumsum(booklets$nit)),
       bk_max_a = booklets$max_a,
       bk_max = lapply(split(design,design$booklet_id), function(d) sum(a[d$last])))
}


# Plausible values

# @param x                summarized respData x, with column pop added if covraite exist
# @param design           tibble(booklet_id <factor>, first <int>, last <int>  ordered by first
# @param b                vector of b's per item_score, including 0 category, ordered by item_id, item_score
#                           or a matrix where each row is the above
# @param a                vector of weights per item_score, including 0 category, ordered by item_id, item_score
# @param nPV              number of plausible values to return per person
# @param from             warm_up: first warm_up pv samples are discarded
# @param step             every step-th sample of pv's is used, typically this would be set to 1 for CML 
#                          item parameters
#                         if NULL, warm_up and step defaults are chosen based on prior and calibration type 
# @param prior_dist       Prior distribution
# @returns                data.frame
#
#


# to do: missing elements in booklet factor or item factor: can they occur and do they mess everything up?
# to do: mixture pv has undesirable ordering with sumscore

pv_chain = function(x, design, b, a, nPV, 
                     warm_up = NULL,
                     step = NULL,
                     A=a, prior_dist = c("normal", "mixture"))
{
  prior_dist = match.arg(prior_dist)
  cal_type = if.else(is.null(dim(b)) || nrow(b)==1,'cml','Bayes')
  pb = get_prog_bar()
  on.exit({pb$close()})

  # klein probleem wordt dat je met deze warmup al snel te weinig iteraties hebt in bayes om 
  # 8 of 16 ketens te voorzien
  # wellicht warm up met 1 bayes param?
  if(prior_dist == 'mixture')
  {
    warm_up = coalesce(warm_up, 150L)
  } else
  {
    warm_up = coalesce(warm_up, if.else(cal_type=='cml', 11L, 20L))
  }
  
  step = coalesce(step, if.else(cal_type=='cml', 1L, 5L))
  
  
  # will make an option
  # let's take 10 chains if we can, but no more than available minus 1, with a minimum of 1 chain
  nchains = min(10L, max(2L, ncores() - 1L))
  
  # for now, we take no more chains than pv's, to not unduly slow the process
  nchains = min(nchains,nPV)
  
  if(cal_type=='Bayes')
  {
    b = t(b)
    iter_per_chain = warm_up + step * (ceiling(nPV/nchains)-1)
    # an option would be to lessen the number of chains
    if(iter_per_chain*nchains > ncol(b)) 
    {
      stop_(sprintf("Need at least %i samples of Bayes NRM parameters to draw %i plausible values on %i cpu-cores", 
                    iter_per_chain*nchains, nPV, nchains))
    }
  } else
  {
    b = matrix(b, ncol=1)
  }
  
  ### prepare bookkeeping
  
  stopifnot(is.factor(x$booklet_id))
  stopifnot(is.factor(design$booklet_id))
  stopifnot(all(levels(x$booklet_id) == levels(design$booklet_id)))
  
  design = pv_design(design, a)

  has_groups = prior_dist=='normal' && 'pop' %in% colnames(x)
  
  if(has_groups)
    x = lapply(split(x, list(x$pop, x$booklet_id), drop = TRUE), function(y) arrange(y, .data$booklet_score))
  else
    x = lapply(split(x, x$booklet_id, drop = TRUE), function(y) arrange(y, .data$booklet_score))

  names(x) = NULL

  
  scoretab = lapply(x, function(y) score_tab_single(y$booklet_score, design$bk_max[[as.character(y$booklet_id[1])]]))
  
  scoretab_counts = tibble(
    booklet_c = sapply(x,function(y) as.integer(y$booklet_id[1])) -1L,
    pop = if.else(has_groups, sapply(x,function(y) y$pop[1]), 1L),
    pop_c = pop-1L,
    n_persons = sapply(x, nrow),
    n_scores = sapply(scoretab, length))
  
  scoretab = unlist(scoretab)

  ### end prepare bookkeeping
  
  # to do: named arguments in cpp call, unnamed is a bit error prone
  
  if(prior_dist=='normal')
  {
    npop = max(scoretab_counts$pop_c) + 1L
    
    start_mu = matrix(rnorm(npop*nchains), ncol=nchains)
    start_sigma = runif(nchains,3,4)
    
    res = pv_chain_normal(b,a,A, 
                          design$first_c, design$last_c, design$bk_cnit, design$bk_max_a,
                          scoretab, scoretab_counts$booklet_c, scoretab_counts$pop_c,
                          scoretab_counts$n_scores, scoretab_counts$n_persons,
                          start_mu, start_sigma, as.integer(nPV), 
                          pb$cpp_prog_init(),
                          as.integer(warm_up), as.integer(step))
    
    dimnames(res$prior_log) = list(var=c('mu','sigma','tau', sprintf("theta_%i",1:(nrow(res$prior_log)-3))),iter=NULL,chain=NULL)
  } else
  {
    start_p = runif(nchains, .4, .6)
    start_mu = matrix(rnorm(2*nchains), nrow=2)
    start_sigma = matrix(runif(2*nchains,1,2), nrow=2)

    res = pv_chain_mix(b,a,A, 
                       design$first_c, design$last_c, design$bk_cnit, design$bk_max_a,
                       scoretab, scoretab_counts$booklet_c, 
                       scoretab_counts$n_scores, scoretab_counts$n_persons,
                       start_mu, start_sigma, start_p, as.integer(nPV), 
                       pb$cpp_prog_init(),
                       as.integer(warm_up), as.integer(step))
    
    dimnames(res$prior_log) = list(var=c('p','mu_1','mu_2','sigma_1','sigma_2'),iter=NULL,chain=NULL)
  }
    
  assign("prior_log", res$prior_log, envir = .GlobalEnv)
  
  colnames(res$theta) = sprintf("PV%i",1:ncol(res$theta))
  
  bind_cols(bind_rows(x), res$theta)
}



# user interface ----------------------------------------------------------


# to do: it only becomes clear in the anon function how many b iterations in bayes we would have needed.
# Either we have to be more lenient in using the same par draw multiple times if necessary
# or we have to shift some of the code to this function

#' Draw plausible values(dev version for testing)
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
#' @param prior_dist use a normal prior for the plausible values or a mixture of two normals. 
#' A mixture is only possible when there are no covariates.
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible values are generated per person (TRUE) or per booklet (FALSE)
#' @param implementation for testing of dev version
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPV plausible values
#' named PV1...PVn.
#' 
#' @details
#' 
#' When the item parameters are estimated using \code{fit_enorm(..., method='Bayes')} and parms_draw = 'sample', the uncertainty 
#' of the item parameters estimates is taken into account when drawing multiple plausible values. 
#' 
#' In there are covariates, the prior distribution is a hierarchical normal with equal variances across groups. When there is only
#' one group this becomes a regular normal distribution. When there are no covariates and prior_dist = "mixture", the prior is a mixture
#' distribution of two normal distributions which gives a little more flexibility than a normal prior. 
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
plausible_values2 = function(dataSrc, parms=NULL, predicate=NULL, covariates=NULL, 
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
  
  plausible_values2_(dataSrc, parms, qtpredicate=qtpredicate, covariates=covariates, nPV=nPV, 
                     parms_draw = parms_draw, env=env,prior_dist = prior_dist ,
                     merge_within_persons=merge_within_persons) %>%
    mutate_if(is.factor, as.character) %>%
    df_format()
}

# to~do: ignore covariate when (some) groups contain to few, <5 say, persons. Add warning.
# what if these are 4 persons with score 0 on an easy test?
# would, in general, the proper way to deal with the pathological case be to add a dummy covariate
# based on characteristics of scoretab? (per booklet and per user covariate of course) 

plausible_values2_ = function(dataSrc, parms=NULL, qtpredicate=NULL, covariates=NULL, nPV=1, parms_draw = c('sample','average'), 
                              env=NULL, prior_dist = c("normal", "mixture"),
                              merge_within_persons=merge_within_persons)
{
  if(is.null(env)) env = caller_env()
  
  if(is.numeric(parms_draw)) parms_draw = as.integer(parms_draw)
  else parms_draw = match.arg(parms_draw)

  prior_dist = match.arg(prior_dist)
  if(!is.null(covariates) && prior_dist == "mixture")
  {
    message('A mixture prior cannot be used together with covariates, setting `prior_dist = "normal"`')
    prior_dist = 'normal'
  }
  

  pb = get_prog_bar(nsteps=if(is.null(parms)) 120 else 100, 
                    retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  if(is.null(parms))
  {
    respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env)
    pb$new_area(20)
    parms = fit_enorm_(respData, method = 'Bayes') 
    
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
  

  if(!is.null(covariates))
  {
    group_number = (function(){i = 0L; function() i <<- i+1L })()
    respData$x = respData$x %>% 
      group_by_at(covariates) %>%
      mutate(pop__ = group_number()) %>%
      ungroup() 
  } 
  
  # join design with the params
  # these can have different levels
  design = suppressWarnings(respData$design %>%
                              inner_join(parms$items, by='item_id') %>% 
                              arrange(.data$booklet_id, .data$first))
  
  
  y = pv_chain(select(respData$x, 'booklet_id', 'person_id', 'booklet_score', any_of(c(pop = 'pop__'))),
               design, parms$b, parms$a, nPV=nPV, prior_dist=prior_dist) %>%
        select(-any_of('pop'))
  

  if(!is.null(covariates))
  {
    # added unique so that booklet_id can be used as a covariate
    y = inner_join(respData$x[,unique(c('booklet_id','person_id',covariates))], y, 
                   by=c('booklet_id','person_id') )
  }
  
  y
}







