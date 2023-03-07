



##########################################
#' Fit the extended nominal response model
#'
#' Fits an Extended NOminal Response Model (ENORM) using conditional maximum likelihood (CML)
#' or a Gibbs sampler for Bayesian estimation.
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param fixed_params Optionally, a prms object from a previous analysis or 
#' a data.frame with parameters, see details.
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nDraws Number of Gibbs samples when estimation method is Bayes. 
#' @param merge_within_persons whether to merge different booklets administered to the same person, enabling linking over persons as well as booklets.
#' @return An object of type \code{prms}. The prms object can be cast to a data.frame of item parameters 
#' using function `coef` or used directly as input for other Dexter functions.
#' @details
#' To support some flexibility in fixing parameters, fixed_params can be a dexter prms object of a data.frame.
#' If a data.frame, it should contain the columns item_id, item_score and a difficulty parameter. Three types of parameters are supported:
#' \describe{
#' \item{delta/beta}{ thresholds between subsequent item categories }
#' \item{eta}{item-category parameters} 
#' \item{b}{exp(-eta)}
#' }
#' Each type corresponds to a different parametrization of the model. 
#' 
#' @references 
#' Maris, G., Bechger, T.M. and San-Martin, E. (2015) A Gibbs sampler for the (extended) marginal Rasch model. 
#' Psychometrika. 80(4), 859-879. 
#' 
#' Koops, J. and Bechger, T.M. and Maris, G. (in press); Bayesian inference for multistage and other 
#' incomplete designs. In Research for Practical Issues and Solutions in Computerized Multistage Testing.
#' Routledge, London. 
#' 
#' @seealso functions that accept a prms object as input: \code{\link{ability}}, \code{\link{plausible_values}}, 
#' \code{\link{plot.prms}}, and \code{\link{plausible_scores}}
#'
fit_enorm = function(dataSrc, predicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), 
                     nDraws=1000, merge_within_persons=FALSE)
{
  method = match.arg(method)
  check_dataSrc(dataSrc)
  check_num(nDraws, 'integer', .length=1, .min=1)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()

  fit_enorm_(dataSrc, qtpredicate = qtpredicate, fixed_params = fixed_params, 
             method=method, nDraws=nDraws, env=env, merge_within_persons=merge_within_persons)
}


fit_enorm_ = function(dataSrc, qtpredicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), 
                      nDraws=1000, env=NULL, merge_within_persons=FALSE) 
{
  method = match.arg(method)
  if(is.null(env)) env = caller_env()

  pb = get_prog_bar(retrieve_data = is_db(dataSrc), nsteps = if.else(method=='Bayes',nDraws,NULL))
  on.exit({pb$close()})

  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env, retain_person_id=FALSE,
                           merge_within_persons = merge_within_persons)
  
  pb$tick()
  ss = get_sufStats_nrm(respData)

  ssI = ss$ssI
  ssIS = ss$ssIS
  design = ss$design
  plt = ss$plt
  scoretab = ss$scoretab
  
  fixed_b = NULL
  has_fixed_parms = !is.null(fixed_params)

  ## deal with fixed parameters
  if(has_fixed_parms)
  {
    if(inherits(fixed_params,'prms'))
    {
      if (fixed_params$inputs$method!="CML")
        message("Posterior means are taken as values for fixed parameters")
      
      fixed_params = fixed_params$inputs$ssIS %>%
        add_column(b = if.else(fixed_params$inputs$method=="CML", fixed_params$est$b, colMeans(fixed_params$est$b)))
      
    } else
    {
      # transform the fixed params to the b parametrization dexter uses internally
      fixed_params = transform.df.parms(fixed_params, out.format = 'b', include.zero = TRUE) 
    }  
    # avoid factor warnings and reduce
    fixed_params$item_id = ffactor(as.character(fixed_params$item_id), levels=levels(design$item_id))
    fixed_params = filter(fixed_params,!is.na(.data$item_id)) 
    
    # check for missing categories in fixed_params, necessary?
    missing_cat = ssIS %>% 
      semi_join(fixed_params, by='item_id') %>%
      left_join(fixed_params, by=c('item_id','item_score')) %>%
      filter(is.na(.data$b) & .data$item_score != 0) 
      
    if(nrow(missing_cat) > 0)
    {
      cat(paste('Some score categories are fixed while some are not, for the same item.',
                'Dexter does not know how to deal with that.\nThe following score categories are missing:\n'))
      missing_cat %>% 
        select(.data$item_id, .data$item_score) %>%
        arrange(.data$item_id, .data$item_score) %>%
        as.data.frame() %>%
        print()
      stop('missing score categories for fixed items, see output')
    }
      
    fixed_b = fixed_params %>%
      right_join(ssIS, by=c('item_id','item_score')) %>%
      arrange(.data$item_id,.data$item_score) %>%
      pull(.data$b)
    
    if(!anyNA(fixed_b)) stop('nothing to calibrate, all parameters are fixed')
  }
  
  if (method=="CML"){
    result = calibrate_CML(scoretab=scoretab, design=design, sufI=ssIS$sufI, a=ssIS$item_score, 
                               first=ssI$first, last=ssI$last, 
                               fixed_b=fixed_b)

  } else 
  {
    result = calibrate_Bayes(scoretab=scoretab, design=design, sufI=ssIS$sufI, a=ssIS$item_score,
                              first=ssI$first, last=ssI$last, nIter=nDraws, fixed_b=fixed_b)
    
    
  }

  mle = design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = theta_MLE(if.else(is.matrix(result$b),colMeans(result$b), result$b), 
                      a=ssIS$item_score, .$first, .$last, se=FALSE)
      theta = est$theta[2:(length(est$theta)-1)]
      tibble(booklet_score=1:length(theta), theta = theta)
    }) %>%
    ungroup() 
  
  outpt = list(est=result, 
               inputs=list(scoretab=scoretab, design=design, ssIS=ssIS, ssI=ssI, design=design,plt=plt,
                            method=method, has_fixed_parms = has_fixed_parms), 
               abl_tables = list(mle = mle),
               xpr=deparse(qtpredicate))
  class(outpt) = append('prms', class(outpt)) 
  outpt
}



# to~do: is there a better plot for calibrate bayes?

#' Plot for the extended nominal Response model
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_enorm
#' @param item_id which item to plot, if NULL, one plot for each item is made
#' @param dataSrc data source, see details
#' @param predicate an expression to subset data in dataSrc
#' @param nbins number of ability groups
#' @param ci confidence interval for the error bars, between 0 and 1. Use 0 to suppress the error bars.
#' Default = 0.95 for a 95\% confidence interval
#' @param ... further arguments to plot
#' @return 
#' Silently, a data.frame with observed and expected values possibly useful to create a numerical fit measure.
#' @details
#' The standard plot shows the fit against the sample on which the parameters were fitted. If
#' dataSrc is provided, the fit is shown against the observed data in dataSrc. This may be useful 
#' for plotting the fit in different subgroups as a visual test for item level DIF. The confidence 
#' intervals denote the uncertainty about the predicted pvalues within the ability groups for the 
#' sample size in dataSrc (if not NULL) or the original data on which the model was fit.
#' 
#' @method plot prms
#' 
plot.prms = function(x, item_id=NULL, dataSrc=NULL, predicate=NULL, nbins=5, ci = .95, ...)
{
  check_num(nbins,'integer',.length=1, .min=2)
  dots = list(...)
  

  if(is.null(item_id))
  {
    if('items' %in% names(dots))
    {
      # common typo
      item_id = dots$items
      dots$items = NULL
    } else
    {
      item_id = x$inputs$ssI$item_id
    }
  }
  if(length(setdiff(item_id,x$inputs$ssI$item_id))>0)
  {
    message('The following items were not found in yourfit object')
    print(setdiff(item_id,x$inputs$ssI$item_id))
    stop('unknown item',call.=FALSE)
  }

  if(!is.null(dataSrc))
  {
    check_dataSrc(dataSrc)
    qtpredicate = eval(substitute(quote(predicate)))
    env = caller_env()
    respData = get_resp_data(dataSrc, qtpredicate, env=env, retain_person_id=FALSE,
                             parms_check=filter(x$inputs$ssIS, .data$item_id %in% local(item_id)))
    
    if(length(setdiff(as.character(item_id), levels(respData$design$item_id)))>0)
    {
      message('The following items were not found in dataSrc')
      print(setdiff(as.character(item_id), levels(x$design$item_id)))
      stop('unknown item',call.=FALSE)
    }
    
    x$abl_tables = list() # in dexmst it was briefly an environment
    
	  x$abl_tables$mle = suppressWarnings({inner_join(respData$design,x$inputs$ssI,by='item_id')}) %>%
      group_by(.data$booklet_id) %>%
      do({
        est = theta_MLE(x$est$b, a=x$inputs$ssIS$item_score, .$first, .$last, se=FALSE)
        theta = est$theta[2:(length(est$theta)-1)]
        tibble(booklet_score=1:length(theta), theta = theta)
      }) %>%
      ungroup() 
    
    x$inputs$plt = get_sufStats_nrm(respData, check_sanity=FALSE)$plt
  }
  

  #many plots
  if(length(item_id) > 1)
  {
    out = lapply(item_id, function(itm) do.call(plot, append(list(x=x, item_id=itm, nbins=nbins, ci=ci), dots)))
    names(out) = as.character(item_id)
    
    return(invisible(out))
  }
  # for dplyr
  item_id_ = as.character(item_id)
  
  expf = expected_score(x, items = item_id)

  max_score = x$inputs$ssIS %>%
    filter(.data$item_id == item_id_) %>%
    pull(.data$item_score) %>%
    max()
  
  plt = x$inputs$plt %>%
    filter(.data$item_id==item_id_) %>%
    inner_join(x$abl_tables$mle, by=c('booklet_id','booklet_score')) %>%
    mutate(abgroup = weighted_ntile(.data$theta, .data$N, n = nbins)) %>%
    group_by(.data$abgroup) %>%
    summarize(gr_theta = weighted.mean(.data$theta,.data$N), avg_score = weighted.mean(.data$meanScore,.data$N), n=sum(.data$N)) %>%
    ungroup() %>%
    mutate(expected_score = expf(.data$gr_theta))
  
  rng = max(plt$gr_theta) - min(plt$gr_theta)
  rng = c(min(plt$gr_theta)-.5*rng/nbins,
          max(plt$gr_theta)+.5*rng/nbins)
  
  plot.args = merge_arglists(dots,
                             default=list(bty='l',xlab = expression(theta), ylab='score',main=item_id),
                             override=list(x = rng,y = c(0,max_score), type="n"))
  
  plot.args$main = fstr(plot.args$main, list(item_id=item_id))
  plot.args$sub = fstr(plot.args$sub, list(item_id=item_id))
  
  do.call(plot, plot.args)
  lines(plt$gr_theta,plt$expected_score, col='grey80') 
  
  plt$outlier = FALSE
  
  if(!is.null(ci) && !is.na(ci) && ci !=0)
  {
    if(ci>1 && ci<100)
      ci = ci/100
    
    if(ci<0 || ci >= 1)
      stop('confidence interval must be between 0 and 1')
    
    qnt = abs(qnorm((1-ci)/2))
    
    I=information(x, items = item_id)
    plt = plt %>%
      mutate(se = sqrt(I(.data$gr_theta)/.data$n),
             conf_min = pmax(.data$expected_score - qnt*.data$se,0),
             conf_max = pmin(.data$expected_score + qnt*.data$se,max_score)) %>%
      mutate(outlier = .data$avg_score < .data$conf_min | .data$avg_score > .data$conf_max)
    
    suppressWarnings({
      arrows(plt$gr_theta, plt$conf_min, 
             plt$gr_theta, plt$conf_max, 
             length=0.05, angle=90, code=3, col='grey80')})
  } 
  
  lines(plt$gr_theta,plt$avg_score)  
  points(plt$gr_theta, plt$avg_score, bg = if_else(plt$outlier, qcolors(1), 'transparent'), pch=21)
  invisible(df_format(plt))
}



print.prms = function(x, ...){
  p = paste0( 'Parameters for the Extended Nominal Response Model\n\n',
              'Method: ', x$inputs$method, ', ',
              ifelse(x$inputs$method == 'CML',
                     paste0('converged in ',x$est$n_iter, ' iterations'),
                     paste0('number of Gibbs samples: ',nrow(x$est$beta))),
              '\nitems: ', nrow(x$inputs$ssI), 
              '\nresponses: ', sum(x$inputs$ssIS$sufI),'\n\n',
              'Use coef() or coefficients() to extract the item parameters.\n')
    
  cat(p)
  invisible(x)
}



#' extract enorm item parameters
#' 
#' @param object an enorm parameters object, generated by the function \code{\link{fit_enorm}}
#' @param hpd width of Bayesian highest posterior density interval around mean_beta, 
#'  value must be between 0 and 1, default is 0.95 
#' @param what which coefficients to return. Defaults to `items` (the item parameters). Can also be `var` for the 
#' variance-covariance matrix (CML only) or `posterior` for all draws of the item parameters (Bayes only)  
#' @param ... further arguments to coef are ignored
#'  
#' @return 
#' Depends on the calibration method and the value of 'what'. For 'items' 
#' 
#' \describe{
#' \item{CML}{a data.frame with columns: item_id, item_score, beta, SE_beta}
#' \item{Bayes}{a data.frame with columns: item_id, item_score, mean_beta, SD_beta, <hpd_b_left>, <hpd_b_right>}
#' }
#' 
#' The posterior distribution and variance covariance matrix are returned as matrices.
#' 
#' 
#' 
coef.prms = function(object, hpd = 0.95, what=c('items','var','posterior'), ...)
{
  x = object
  what = match.arg(what)
  
  if(what=='items')
  {
    if (x$inputs$method=="CML")
    {
      atab=data.frame(item_id=x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                      item_score=as.integer(x$inputs$ssIS$item_score[-x$inputs$ssI$first]),
                      beta=x$est$beta,
                      SE_beta=sqrt(diag(x$est$acov.beta)),stringsAsFactors=FALSE)
    } else
    {
      if(hpd <= 0 ||  hpd >= 1)
        stop('hpd must be between 0 and 1')
      
      hh = t(apply(x$est$beta,2,hpdens, conf=hpd))
      atab=data.frame(item_id = x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                      a = x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                      mb = colMeans(x$est$beta),
                      sdb = apply(x$est$beta, 2, sd),
                      hpdl = hh[,1], hpdr=hh[,2],stringsAsFactors=FALSE)
      colnames(atab)=c("item_id" ,"item_score", "mean_beta", "SD_beta", 
                       sprintf("%i_hpd_b_left", round(100 * hpd)),
                       sprintf("%i_hpd_b_right", round(100 * hpd)))
    }
    atab$item_id = as.character(atab$item_id)
    rownames(atab) = NULL
    return(df_format(atab))
  } else if(what=='var')
  {
    if(x$inputs$method!="CML")
      stop('Variance-covariance matrix is only available for CML extimation')
    m = x$est$acov.beta
    colnames(m) = rownames(m) = paste(x$inputs$ssIS$item_id, x$inputs$ssIS$item_score)[-x$inputs$ssI$first]
    return(m)
  } else if(what=='posterior')
  {
    if(x$inputs$method!="Bayes")
      stop('The posterior of item parameters is only available for Bayesian estimation')
    m=x$est$beta
    colnames(m) = paste(x$inputs$ssIS$item_id, x$inputs$ssIS$item_score)[-x$inputs$ssI$first]
    return(m)
  }
}


#' Functions of theta
#' 
#' returns information function, expected score function, score simulation function, or score distribution 
#' for a single item, an arbitrary group of items or all items
#' 
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and, 
#' depending on parametrization, a column named either beta/delta, eta or b
#' @param items vector of one or more item_id's. If NULL and booklet_id is also NULL, all items in parms are used
#' @param booklet_id id of a single booklet (e.g. the test information function), if items is not NULL this is ignored
#' @param which.draw the number of the random draw (only applicable if calibration method was Bayes). If NULL, the mean 
#' beta parameter will be used
#' 
#' @return Each function returns a new function which accepts a vector of theta's. These return the following values: 
#' \describe{
#' \item{information}{an equal length vector with the information estimate at each value of theta.}
#' \item{expected_score}{an equal length vector with the expected score at each value of theta}
#' \item{r_score}{a matrix with length(theta) rows and one column for each item containing simulated scores based on theta. 
#' To obtain test scores, use rowSums on this matrix}
#' \item{p_score}{a matrix with length(theta) rows and one column for each possible sumscore containing the probability of 
#' the score given theta}
#' }
#' 
#' @examples
#' 
#' db = start_new_project(verbAggrRules,':memory:')
#' add_booklet(db,verbAggrData, "agg")
#' p = fit_enorm(db)
#' 
#' # plot information function for single item
#' 
#' ifun = information(p, "S1DoScold")
#' 
#' plot(ifun,from=-4,to=4)
#' 
#' # compare test information function to the population ability distribution
#' 
#' ifun = information(p, booklet="agg")
#' 
#' pv = plausible_values(db,p)
#' 
#' op = par(no.readonly=TRUE)
#' par(mar = c(5,4,2,4))
#' 
#' plot(ifun,from=-4,to=4, xlab='theta', ylab='test information')
#' 
#' par(new=TRUE)
#' 
#' plot(density(pv$PV1), col='green', axes=FALSE, xlab=NA, ylab=NA, main=NA)
#' axis(side=4)
#' mtext(side = 4, line = 2.5, 'population density (green)')
#' 
#' par(op)
#' close_project(db)
#' 
information = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='information')
}

#' @rdname information
expected_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='expected')
}

#' @rdname information
r_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='sim')
}

#' @rdname information
p_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='pmf')
}


theta_function = function(parms, items=NULL, booklet=NULL, which.draw=NULL, 
                          what=c('information','expected','sim','pmf'))
{
  what = match.arg(what)
  check_character(items,nullable=TRUE)
  check_string(booklet,name='booklet_id',nullable=TRUE)
  check_num(which.draw,nullable=TRUE)
  
  # data preparation
  # create fl(item_id,first,last), a, b
  
  if(inherits(parms,'data.frame'))
  {
    out = transform.df.parms(parms,'b',include.zero=TRUE)
    a = out$item_score
    b = out$b
    
    fl = out %>%
      mutate(rn=row_number()) %>%
      group_by(.data$item_id) %>%
      summarize(first=as.integer(min(.data$rn)), last=as.integer(max(.data$rn))) %>%
      ungroup()
    
    if(!is.null(items))
      fl = filter(fl, fl$item_id %in% items)
    
  } else if(inherits(parms,'prms'))
  {
    a = parms$inputs$ssIS$item_score
    b = parms$est$b
    if(is.matrix(b))
    {
      if(is.null(which.draw))
      {
        b = colMeans(b)
      } else
      {
        if(which.draw<1 || which.draw > nrow(b))
          stop('argument `which.draw` out of range')
        b = as.vector(b[which.draw,])
      }
    }
    
    fl = parms$inputs$ssI
    fl$item_id = as.character(fl$item_id)
    
    if(!is.null(items))
    {
      items = unique(items)
      suppressWarnings({fl = semi_join(fl, tibble(item_id=items),by='item_id')})
      if(nrow(fl) != length(items))
      {
        message('unknown items:')
        print(sort(setdiff(items,fl$item_id)))
        stop('Some items were not found, see output')
      }
    } else if(!is.null(booklet))
    {
      booklet = unique(booklet)
      design = parms$inputs$design
      if(length(intersect(booklet,design$booklet_id))<length(booklet))
      {
        stop('unknown booklet')
      }
      
      fl = design %>%
        filter(.data$booklet_id %in% booklet) %>%
        distinct(.data$item_id, .keep_all=TRUE)
    }  
    fl = arrange(fl,.data$first)
  } else stop_('parms must be a data.frame or an object of type `prms`') 
  rm(parms)  
  #output
  
  if(what=='information')
  {
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta)))
        stop('theta may not contain nan/NA values')  
      
      res = rep(0,length(theta))
      res[is.finite(theta)] = 
        IJ_(b,a,fl$first, fl$last, theta[is.finite(theta)])$I
      # extremely large values overflow to NaN, recover as 0
      res[is.nan(res)] = 0
      res
    }
    class(out) = append('inf_func',class(out))
    
  } else if(what == 'expected')
  {
    max_score = sum(a[fl$last])
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta))) 
        stop('theta may not contain nan/NA values') 
      
      res = rep(0,length(theta))
      
      res[is.finite(theta)] = E_score(theta[is.finite(theta)],  
                                      b=b, a=a, 
                                      first=fl$first, last=fl$last)
      
      res[is.infinite(theta) & theta > 0] = max_score
      # extremely large values of theta overflow to NaN (small values undeflow to zero, which is fine)
      res[is.nan(res)] = max_score
      res
    }
    class(out) = append('exp_func',class(out))
  } else if(what=='sim')
  {
    out = function(theta)
    {
      res = rscore_item(theta,b=b,a=a,first = fl$first, last = fl$last)
      colnames(res) = fl$item_id
      res
    }
    class(out) = append('sim_func',class(out))
  } else if(what=='pmf')
  {
    out = function(theta)
    {
      res = pscore(theta,b=b,a=a,first = fl$first, last = fl$last)
      t(res)
    }
    class(out) = append('pmf_func',class(out))
  }
  
  out
}

print.inf_func = function(x,...) cat('Information function: I(theta)\n')
print.exp_func = function(x,...) cat('Conditional expected score function: E(X_i|theta)\n')
print.sim_func = function(x,...) cat('function to simulate item scores: (x_i1, ..., x_ip) ~ ENORM(theta)\n')
print.pmf_func = function(x,...) cat('Conditional score distribution function: P(x_+|theta)\n')


#' Latent correlations
#'
#' Estimates correlations between latent traits. Use an item_property to distinguish the different scales. 
#' This function uses plausible values so results may differ slightly between calls. 
#' Note: this is a new and slightly experimental function and therefore still a bit slow. It will probably 
#' become faster in future versions of dexter.
#'
#' @param dataSrc A connection to a dexter database or a data.frame with columns: person_id, item_id, item_score and 
#' the item_property
#' @param item_property The name of the item property used to define the domains. If \code{dataSrc} is a dexter db then the
#' item_property must match a known item property. If datasrc is a data.frame, item_property must be equal to
#'  one of its column names.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param nDraws Number of draws for plausible values
#' @param hpd width of Bayesian highest posterior density interval around the correlations, 
#'  value must be between 0 and 1.
#' @param use Only complete.obs at this time. Respondents who don't have a score for one or more scales are removed.
#' 
#' @return List containing a estimated correlation matrix, the corresponding standard deviations, 
#' and the lower and upper limits of the highest posterior density interval
#' 
latent_cor = function(dataSrc, item_property, predicate=NULL, nDraws=500, hpd=0.95, use="complete.obs")
{
  check_dataSrc(dataSrc)
  check_num(nDraws, 'integer', .length=1, .min=1)
  check_num(hpd,  .length=1, .min=1/nDraws,.max=1)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  # use = pmatch(use, c("all.obs", "complete.obs", 
  #                     "pairwise.complete.obs", "everything", "na.or.complete"))
  
  if(use != "complete.obs")
    stop("only 'complete.obs' is currently implemented")
  
  pb = get_prog_bar(nDraws, retrieve_data = is_db(dataSrc), lock=TRUE)
  on.exit({pb$close()})
  
  from = 5
  by = 2
  which.keep = seq(from,(from-by)*(from>by)+by*nDraws,by=by)
  nIter=max(which.keep)
  
  if(is.matrix(dataSrc))
  {
    # to do: we can take a char vector of length ncol, but that goes for all
    # analyses functions with item_properties, ip's can be made more flexible for df as well
    stop("a matrix datasrc is not yet implemented for this function")
  }
  # to do: this is a bit tricky, we will often need to merge over persons, e.g. if booklets are administered
  # per subject. But if the same person makes multiple tests for the same subject, this is not good.
  # Still have to decide on a solution
  
  respData = get_resp_data(dataSrc, qtpredicate, env=env, extra_columns=item_property,
                           merge_within_persons=TRUE)

  respData$x[[item_property]] = ffactor(as.character(respData$x[[item_property]]))
  lvl = levels(respData$x[[item_property]])
  nd = length(lvl)
  
  pb$set_nsteps(nIter + 4*nd)
  
  respData$x = respData$x %>%
    group_by(.data$person_id) %>%
    filter(nd == n_distinct(.data[[item_property]])) %>%
    ungroup()
  
  respData$x$person_id = ffactor(respData$x$person_id,as_int=TRUE)
  
  
  np = max(respData$x$person_id)
  respData = lapply(split(respData$x, respData$x[[item_property]]), get_resp_data)
  models = lapply(respData, function(x){try(fit_enorm(x), silent=TRUE)})
  names(models)=names(respData)
  if(any(sapply(models,inherits,what='try-error')))
  {
    message('\nThe model could not be estimated for one or more item properties, reasons:')
    models[!sapply(models,inherits,what='try-error')] = 'OK'
    print(data.frame(item_property=names(models), 
                      result=gsub('^.+try\\(\\{ *:','',trimws(as.character(models)))))
    stop('Some models could not be estimated')
  }
  
  abl = mapply(ability, respData, models, SIMPLIFY=FALSE,
               MoreArgs = list(method="EAP", prior="Jeffreys",standard_errors=FALSE))
  
  pb$tick(2*nd)
  
  # some matrices
  # we cannot rely on order in ability so this is the way for now
  abl_mat = matrix(NA_real_,np,nd)
  for(d in 1:nd)
    abl_mat[abl[[d]]$person_id,d] = abl[[d]]$theta

  
  acor = cor(abl_mat,use=use)
  pv = matrix(0,np,nd)
  
  reliab = rep(0,nd)
  sd_pv = rep(0,nd)
  mean_pv = rep(0,nd)
  for (i in 1:nd)
  {
    pvs = plausible_values(respData[[i]],models[[i]],nPV = 2)
    reliab[i] = cor(pvs$PV1,pvs$PV2)
    sd_pv[i] = sd(pvs$PV1)
    mean_pv[i] = mean(pvs$PV1)
    pv[pvs$person_id,i] = pvs$PV1
  }
  pb$tick(2*nd)
  
  for (i in 1:(nd-1))
  {
    for (j in ((i+1):nd)){
      acor[i,j] = acor[i,j]/sqrt(reliab[i]*reliab[j])
      acor[j,i] = acor[i,j]
    }
  }
  
  out_sd=matrix(0,nd,nd)
  out_cor = acor
  
  # make everything simple and zero indexed
  models = lapply(models, simplify_parms, zero_indexed=TRUE)
  respData = lapply(respData, get_resp_data, summarised=TRUE)
  for(d in 1:nd)
  {
    respData[[d]]$x$booklet_id = as.integer(respData[[d]]$x$booklet_id) - 1L
    models[[d]]$bcni = c(0L,cumsum(table(as.integer(models[[d]]$design$booklet_id))))
  }
  
  prior = list(mu=mean_pv, Sigma = diag(1.7*sd_pv) %*% acor %*% diag(1.7*sd_pv))
  
  store = matrix(0, length(which.keep), length(as.vector(prior$Sigma)))
  tel = 1
  for (i in 1:nIter)
  {
    for (d in 1:nd)
    {
      cons = condMoments(prior$mu, prior$Sigma, d, x.value=pv[,-d]) 

      PV_sve(models[[d]]$b, models[[d]]$a, models[[d]]$design$first, models[[d]]$design$last, 					
             models[[d]]$bcni,
             respData[[d]]$x$booklet_id, respData[[d]]$x$booklet_score, cons$mu, sqrt(cons$sigma),
             pv,d-1L, 10L)
    }
    
    
    prior = update_MVNprior(pv,prior$Sigma)

    if (i %in% which.keep){
      store[tel,] = as.vector(cov2cor(prior$Sigma))
      tel=tel+1
    }
    pb$tick()
  }
  out_sd = matrix(apply(store,2,sd),nd,nd)
  diag(out_sd)=0
  pd = t(apply(store,2,hpdens, conf=hpd))
  res = list(cor = matrix(colMeans(store),nd,nd), sd = out_sd,
             hpd_l = matrix(pd[,1], nd, nd), hpd_h = matrix(pd[,2], nd, nd))

  lapply(res, function(x){colnames(x) = rownames(x) = names(models); x})
}
