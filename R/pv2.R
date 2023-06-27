
eat_vector = function(x)
{
  i=1L
  function(n=1L)
  {
    start = i
    i <<- n + i
    x[start:(i-1L)]
  }
}

# combine some bookkeeping for pv_normal and pv_mixture
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

# @param x                summarized respData x with column pop added
# @param design           tibble(booklet_id <factor>, first <int>, last <int>  ordered by first
# @param b                vector of b's per item_score, including 0 category, ordered by item_id, item_score
#                           or a matrix where each row is the above
# @param a                vector of weights per item_score, including 0 category, ordered by item_id, item_score
# @param nPV              number of plausible values to return per person
# @param from             warm_up: first warm_up pv samples are discarded
# @param step             every step-th sample of pv's is used, typically this would be set to 1 for CML 
#                          item parameters
# @param prior_dist       Prior distribution
# @returns
#
#

pv_normal = function(x, design, b, a, nPV, 
               cal_type = ifelse(is.null(dim(b)) || nrow(b)==1,'cml','Bayes'),
               warm_up = ifelse(cal_type=='Bayes', Gibbs.settings$from.pv.bayes, Gibbs.settings$from.pv.cml),
               step = ifelse(cal_type=='Bayes', Gibbs.settings$step.pv.bayes, Gibbs.settings$step.pv.cml),
               A=a)
{
  cal_type = match.arg(cal_type, c('cml','Bayes'))

  ### prepare iterations
  
  # matrix to fill
  theta = matrix(0, nrow(x), nPV)
  colnames(theta) = sprintf('PV%i',1:nPV)
  
  which.pv = seq(warm_up,(warm_up-step)*(warm_up > step) + step*nPV, by=step)
  niter = max(which.pv)

  if(cal_type=='Bayes')
  {
    b = t(b)
    bstep = as.integer((ncol(b)-1L)/(niter-1L))
  } else
  {
    b = matrix(b, ncol=1)
    bstep = 0L
  }

  
  ### end prepare iterations
  
  
  ### prepare bookkeeping
  
  stopifnot(is.factor(x$booklet_id))
  stopifnot(is.factor(design$booklet_id))
  stopifnot(all(levels(x$booklet_id) == levels(design$booklet_id)))

  design = pv_design(design, a)
  #arrange person_id is for testing
  x = lapply(split(x, list(x$pop, x$booklet_id), drop = TRUE), function(y) arrange(y, .data$booklet_score, person_id))
  names(x) = NULL
  pop_indx = unlist(lapply(x,'[[', 'pop'))
  
  scoretab = lapply(x, function(y) score_tab_single(y$booklet_score, design$bk_max[[as.character(y$booklet_id[1])]]))
  
  scoretab_counts = tibble(
    booklet_c = sapply(x,function(y) as.integer(y$booklet_id[1])) -1L,
    pop = sapply(x,function(y) y$pop[1]),
    pop_c = pop-1L,
    n_persons = sapply(x, nrow),
    n_scores = sapply(scoretab, length))
  
  scoretab = unlist(scoretab)
  cscoretab = cumsum(scoretab)
  
  npop = max(scoretab_counts$pop_c) + 1L
  
  # cumulative tab scores for loop
  cn_scores = c(0L,cumsum(scoretab_counts$n_scores))
  
  ### end prepare bookkeeping
  
  ### starting values

  priors = list(mu=rep(0,npop), sigma=rep(4,npop), mu.a=0, sigma.a=1)

  ### end starting values

  
  current_pv_col = 1L
  current_b_col = 1L # to get equal to original: current_b_col = bstep

  for(iter in 1:niter)
  {
    pv_draw(b[,current_b_col], a, A,
            design$first_c, design$last_c, design$bk_cnit, design$bk_max_a,
            scoretab, cscoretab, scoretab_counts$booklet_c, scoretab_counts$pop_c,
            scoretab_counts$n_scores,cn_scores, scoretab_counts$n_persons,
            priors$mu, priors$sigma,
            theta, current_pv_col-1L)
    

    if(iter == niter) break
    
    # update prior
    if(npop == 1)
      priors = update_pv_prior(theta[,current_pv_col], pop_indx, priors$mu, priors$sigma)
    else
      priors = update_pv_prior_H(theta[,current_pv_col], pop_indx, priors$mu, priors$sigma, priors$mu.a, priors$sigma.a)
    
    current_b_col = current_b_col + bstep
    if(iter %in% which.pv) current_pv_col = current_pv_col + 1L
  }

  bind_cols(bind_rows(x), theta)
}




pv_mixture = function(x, design, b, a, nPV, 
                      cal_type = ifelse(is.null(dim(b)) || nrow(b)==1,'cml','bayes'),
                      warm_up = ifelse(cal_type=='bayes', Gibbs.settings$from.pv.bayes, Gibbs.settings$from.pv.cml),
                      step = ifelse(cal_type=='bayes', Gibbs.settings$step.pv.bayes, Gibbs.settings$step.pv.cml),
                      A=a)
{
  # still a bit too much copy paste from the previous
  
  cal_type = match.arg(cal_type, c('cml','Bayes'))
  
  ### prepare iterations
  
  # matrix to fill
  theta = matrix(0, nrow(x), nPV)
  colnames(theta) = sprintf('PV%i',1:nPV)
  
  which.pv = seq(warm_up,(warm_up-step)*(warm_up > step) + step*nPV, by=step)
  niter = max(which.pv)
  
  if(cal_type=='Bayes')
  {
    b = t(b)
    bstep = as.integer((ncol(b)-1L)/(niter-1L))
  } else
  {
    b = matrix(b, ncol=1)
    bstep = 0L
  }
  
  ### end prepare iterations
  
  
  ### prepare bookkeeping
  
  stopifnot(is.factor(x$booklet_id))
  stopifnot(is.factor(design$booklet_id))
  stopifnot(all(levels(x$booklet_id) == levels(design$booklet_id)))
  
  design = pv_design(design, a)
  
  x = lapply(split(x, x$booklet_id), function(y) arrange(y, .data$booklet_score))
  names(x) = NULL

  
  make_scoretab = function(x, bkmax, group)
  {
    p = eat_vector(group)
    
    scoretab = lapply(x, 
      function(y) score_tab_double(y$booklet_score, p(nrow(y)), design$bk_max[[as.character(y$booklet_id[1])]]))
    
    scoretab_counts = tibble(
      booklet_c = rep(sapply(x,function(y) as.integer(y$booklet_id[1])) -1L, each=2),
      group = rep(0:1,times = length(x)),
      n_persons = as.integer(sapply(scoretab, colSums)),
      n_scores = rep(sapply(scoretab, nrow), each=2))
    
    scoretab = unlist(lapply(scoretab, as.integer))
    
    
    list(scoretab = scoretab, cscoretab = cumsum(scoretab),
         booklet_c = scoretab_counts$booklet_c, group = scoretab_counts$group,
         n_scores = scoretab_counts$n_scores, cn_scores = c(0L,cumsum(scoretab_counts$n_scores)),
         n_persons = scoretab_counts$n_persons)
  }
  
  priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2), 
                group =  sample(0:1, nrow(theta), replace=TRUE, prob=c(0.5,0.5)))

  ### end prepare bookkeeping
  
  
  current_pv_col = 1L
  current_b_col = 1L

  for(iter in 1:niter)
  {
    stab = make_scoretab(x, bkmax, priors$group)
    
    pv_draw(b[,current_b_col], a, A,
            design$first_c, design$last_c, design$bk_cnit, design$bk_max_a,
            stab$scoretab, stab$cscoretab, stab$booklet_c, stab$group, 
            stab$n_scores,stab$cn_scores, stab$n_persons,
            priors$mu, priors$sigma, 
            theta,current_pv_col-1L) 
    
    theta[,current_pv_col] = arrange_pv(theta[,current_pv_col], stab$n_persons, priors$group)
    
    if(iter == niter) break
    
    # update prior
    
    priors = update_pv_prior_mixnorm(theta[,current_pv_col], priors$p, priors$mu, priors$sigma)
    priors$group = priors$grp -1L
    
    current_b_col = current_b_col + bstep
    if(iter %in% which.pv) current_pv_col = current_pv_col + 1L
  }

  bind_cols(bind_rows(x), theta)
}








# user interface ----------------------------------------------------------


#' Draw plausible values(dev version fro testing)
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
  } else if(prior_dist == 'normal')
  {
    # niet varierende pop toevoegen maakt code in pv eenvoudiger
    respData$x$pop__ = 1L
  }
  
  # join design with the params
  # these can have different levels
  design = suppressWarnings(respData$design %>%
                              inner_join(parms$items, by='item_id') %>% 
                              arrange(.data$booklet_id, .data$first))
  
  if(prior_dist == 'mixture')
  {
    y = pv_mixture(select(respData$x, 'booklet_id', 'person_id', 'booklet_score'),
         design, parms$b, parms$a, nPV=nPV) 
  } else
  {
    y = pv_normal(select(respData$x, 'booklet_id', 'person_id', 'booklet_score', pop = 'pop__'),
            design, parms$b, parms$a, nPV=nPV) %>%
      select(-'pop')
  }

  
  

  if(!is.null(covariates))
  {
    # added unique so that booklet_id can be used as a covariate
    y = inner_join(respData$x[,unique(c('booklet_id','person_id',covariates))], y, 
                   by=c('booklet_id','person_id') )
  }
  
  y
}





