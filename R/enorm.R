
# to do: reparametrization, delta in mirt for polytomous is different

# mpar = lapply(cmirt, function(v)
# {
#   tibble(delta = -v[grepl("^d[^0]*$", colnames(v),perl=T)], item_score=1:length(delta))
# }) %>%
#   bind_rows(.id='item_id') %>%
#   group_by(item_id) %>%
#   mutate(delta = 2* delta - cumsum(delta)) %>%
#   ungroup()


# item_id, item_score, and b, eta of beta
transform.df.parms = function(parms.df, out.format = c('b','beta','eta'), include.zero = TRUE)
{
  # start with many checks
  out.format = match.arg(out.format)
  colnames(parms.df) = tolower(colnames(parms.df))
  
  if('delta' %in% colnames(parms.df))
    parms.df = rename(parms.df, beta = 'delta')
  in.format = intersect(colnames(parms.df), c('b','beta','eta'))
  
  if(length(in.format) == 0)
    stop('parameters must contain at least one of following columns: b, beta, eta')
  
  if(length(in.format)>1)
  {
    in.format = in.format[1]
    message(paste0("Using '",in.format,"' as input parameter"))
  }
  
  if(!all(c('item_id','item_score') %in% colnames(parms.df)))
    stop('parameters must contain the columns: item_id, item_score')
  
  if(any(parms.df$item_score%%1 > 0))
    stop("column 'item_score' must be integer valued")
  
  if(n_distinct(parms.df$item_id, parms.df$item_score) < nrow(parms.df))
    stop('multiple parameters supplied for the same item and score')
  
  parms.df = parms.df %>% 
    mutate(item_id = as.character(.data$item_id), item_score = as.integer(.data$item_score)) %>%
    arrange(.data$item_id, .data$item_score)
  
  
  mm = parms.df %>% 
    group_by(.data$item_id) %>% 
    summarise(min_score = min(.data$item_score), max_score = max(.data$item_score)) %>%
    ungroup()
  
  in.zero = any(mm$min_score == 0)
  
  if(in.zero && any(mm$min_score) != 0)
    stop("Either all items or none of the items should include a zero score parameter")
  
  if(any(mm$max_score == 0))
    stop('All items should contain at least one non-zero score parameter')
  
  if(any(mm$min_score<0))
    stop("Negative scores are not allowed")

  if(in.format == 'b' && any(parms.df$b < 0))
      stop("A 'b' parameter cannot be negative, perhaps you meant to include a 'beta' parameter?")
  
  fl = parms.df %>%
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>% 
    summarize(first = as.integer(min(.data$rn)), last = as.integer(max(.data$rn))) %>%
    ungroup()
  
  args = list(first = fl$first, last = fl$last, parms.df = parms.df, 
              out.zero = include.zero, in.zero = in.zero)
  do.call(get(paste0(in.format,'2',out.format)), args)
}


# to do: print first 30 items without score variation, extremely annoying to get just a message

##########################################
#' Fit the extended nominal response model
#'
#' Fits an Extended NOminal Response Model (ENORM) using conditional maximum likelihood (CML)
#' or a Gibbs sampler for Bayesian estimation.
#'
#'
#' @param dataSrc Data source: a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param fixed_params Optionally, a prms object from a previous analysis or 
#' a data.frame with columns: item_id, item_score (omitting 0 score category) and beta. To facilitate the user in entering parameter values, we assume the parameterisation used by OPLM; in short, beta's are thresholds between categories. At this moment, it is not possible to fix some but not all categories of an item.
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nIterations Number of Gibbs samples when estimation method is Bayes. The maximum 
#' number of iterations when using CML.
#' @param link_persons whether to merge different booklets administered to the same person, enabling linking over persons instead of just booklets.
#' @return An object of type \code{prms}. The prms object can be cast to a data.frame of item parameters 
#' using function `coef` or used directly as input for other Dexter functions.
#' 
#' @references 
#' Maris, G., Bechger, T.M. and San-Martin, E. (2015) A Gibbs sampler for the (extended) marginal Rasch model. 
#' Psychometrika. 2015; 80(4): 859–879. 
#' 
#' @seealso functions that accept a prms object as input: \code{\link{ability}}, \code{\link{plausible_values}}, 
#' \code{\link{plot.prms}}
#'
fit_enorm = function(dataSrc, predicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), 
                     nIterations=500, link_persons=FALSE)
{
  method = match.arg(method)
  check_dataSrc(dataSrc)
  check_num(nIterations, 'integer', .length=1, .min=1)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  fit_enorm_(dataSrc, qtpredicate = qtpredicate, fixed_params = fixed_params, 
             method=method, nIterations=nIterations, env=env, link_persons=link_persons)
}


fit_enorm_ = function(dataSrc, qtpredicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), 
                      nIterations=500, env=NULL, link_persons=FALSE) 
{
  method = match.arg(method)
  if(is.null(env)) env = caller_env()
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env, retain_person_id=FALSE,
                           merge_within_person = link_persons)
  
  design = respData$design
  if(nrow(respData$x) == 0) stop('no data to analyse')
  
  # to do: I think we need tests with at least three items, fischer criterium
  if(nrow(design) == 1) 
    stop('There are responses to only one item in your selection, this cannot be calibrated.') 
  
  if(!is_connected(design))
    stop('Your design is not connected')  

  ss = get_sufStats_nrm(respData)
  ssIS = ss$ssIS
  plt = ss$plt
  
  # to do: this should be tested in get_resp_data, it can really mess up the c code
  if(min(ssIS$item_score)<0)
    stop("item_scores must be positive numbers")
  
  # bug in dplyr, min/max of integer in group_by becomes double
  # to do: check where min/max is used elsewhere
  ssI  = ssIS %>% 
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>%
    summarise(first = as.integer(min(.data$rn)),last = as.integer(max(.data$rn))) %>%
    ungroup() %>%
    arrange(.data$item_id)
  
  if(any(ssI$first == ssI$last)) 
    stop('One or more items are without score variation') # to do, show which
  
  design = design %>%
    inner_join(ssI,by='item_id') %>%
    arrange(.data$booklet_id, .data$first)
  
  itm_max = ssIS %>% 
    group_by(.data$item_id) %>% 
    summarise(maxScore = as.integer(max(.data$item_score))) %>% 
    ungroup()
  
  # max booklet scores
  maxScores = itm_max %>%
    inner_join(design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    summarise(maxTotScore = sum(.data$maxScore))
  
  # booklets 0:maxscore
  all_scores = maxScores %>% 
    group_by(.data$booklet_id) %>%
    do({tibble(booklet_score=0:.$maxTotScore)}) %>%
    ungroup()
  
# to do: check if this is efficient for e.g. stex
  scoretab = plt %>%
    distinct(.data$booklet_id, .data$booklet_score,.data$N) %>%
    right_join(all_scores, by=c('booklet_id','booklet_score')) %>%
    mutate(N=coalesce(.data$N, 0L)) %>%
    arrange(.data$booklet_id, .data$booklet_score)

  #to do: temp check, when this fails it is likely to be the merge over person problem (see stex data)
  if(nrow(all_scores)!=nrow(scoretab))
    stop("mistake in scoretab")
  
  fixed_b = NULL
  has_fixed_parms = !is.null(fixed_params)
  # to do: factors and bayes: warnings?
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
    fixed_params$item_id = ffactor(fixed_params$item_id, levels=levels(design$item_id))
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
    
    # to do: this test fails for some reason
    if(!anyNA(fixed_b)) stop('nothing to calibrate, all parameters are fixed')
  }
  
  if (method=="CML"){
    result = calibrate_CML(scoretab=scoretab, design=design, sufI=ssIS$sufI, a=ssIS$item_score, 
                               first=ssI$first, last=ssI$last, nIter=nIterations,
                               fixed_b=fixed_b)
  } else 
  {
    result = calibrate_Bayes(scoretab=scoretab, design=design, sufI=ssIS$sufI, a=ssIS$item_score,
                              b=exp(runif(nrow(ssIS), -1, 1)), 
                              first=ssI$first, last=ssI$last, nIter=nIterations, fixed_b=fixed_b)
  }
  
  mle = design %>% 
    group_by(.data$booklet_id) %>%
    do({
      est = theta_MLE(result$b, a=ssIS$item_score, .$first, .$last, se=FALSE)
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



# to do: is there a better plot for calibrate bayes?
#' Plot for the extended nominal Response model
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_enorm
#' @param item_id which item to plot, if NULL, one plot for each item is made
#' @param nbins number of ability groups
#' @param ci confidence interval for the error bars, between 0 and 1. 0 means no error bars.
#' Default = 0.95 for a 95\% confidence interval
#' @param ... further arguments to plot
#' 
#' @method plot prms
#' 
# to do: make item_id everywhere, r's partial match will take care of `item`
plot.prms = function(x,item_id=NULL, nbins=5, ci = .95, ...)
{
  if(is.null(x$inputs$plt))
  {
    if(inherits(x,'mst_enorm')) stop('Sorry, the plot method for enorm_mst will only be supported in dexterMST from version 0.1.1 onwards')
    stop('Sorry, the plot method is only available for parameter objects produced with dexter 0.8.1 or later')
  }
  check_num(nbins,'integer',.length=1, .min=2)
  
  # if no item id provided plot them all
  if(is.null(item_id))
    item_id = x$inputs$ssI$item_id

  # this was implemented back when ability_tables mle was kind of slow
  # to do: can be removed in time (~jan 2020), plt can also be pre-processed to include theta
  # remove infinite, et.c
  if(is.null(x$abl_tables$mle))
    x$abl_tables$mle = ability_tables(x, standard_errors=FALSE, method='MLE')
  
  if(length(item_id) > 1)
  {
    return(invisible(
      lapply(item_id, function(itm) plot(x, itm, nbins=nbins, ci=ci, ...))))
  }
  # for dplyr
  item_id_ = item_id
  
  expf = expected_score(x, items = item_id)

  max_score = x$inputs$ssIS %>%
    filter(.data$item_id == item_id_) %>%
    pull(.data$item_score) %>%
    max()
  
  plt = x$inputs$plt %>%
    filter(.data$item_id==item_id_) %>%
    inner_join(x$abl_tables$mle, by=c('booklet_id','booklet_score')) %>%
    filter(is.finite(.data$theta)) %>%
    mutate(abgroup = weighted_ntile(.data$theta, .data$N, n = nbins)) %>%
    group_by(.data$abgroup) %>%
    summarize(gr_theta = weighted.mean(.data$theta,.data$N), avg_score = weighted.mean(.data$meanScore,.data$N), n=sum(.data$N)) %>%
    ungroup() %>%
    mutate(expected_score = expf(.data$gr_theta))
  
  rng = max(plt$gr_theta) - min(plt$gr_theta)
  rng = c(min(plt$gr_theta)-.5*rng/nbins,
          max(plt$gr_theta)+.5*rng/nbins)
  
  plot.args = merge_arglists(list(...),
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
    
    cmin = function(p, n) pmax(0, p - qnt * sqrt(p*(1-p)/n))
    cmax = function(p, n) pmin(1, p + qnt * sqrt(p*(1-p)/n))
    
    plt = plt %>%
      mutate(conf_min = max_score * cmin(.data$expected_score/max_score, .data$n),
             conf_max = max_score * cmax(.data$expected_score/max_score, .data$n)) %>%
      mutate(outlier = .data$avg_score < .data$conf_min | .data$avg_score > .data$conf_max)
    
    arrows(plt$gr_theta, plt$conf_min, 
           plt$gr_theta, plt$conf_max, 
           length=0.05, angle=90, code=3, col='grey80')
  } 
  
  lines(plt$gr_theta,plt$avg_score)  
  points(plt$gr_theta, plt$avg_score, bg = if_else(plt$outlier, qcolors(1), 'transparent'), pch=21)
  invisible(plt)
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
#' @param bayes_hpd_b width of Bayesian highest posterior density interval around mean_beta, 
#'  value must be between 0 and 1, default is 0.95 
#' @param ... further arguments to coef are ignored
#'  
#' @return 
#' Depends on the calibration method:
#' \describe{
#' \item{for CML}{a data.frame with columns: item_id, item_score, beta, SE_beta}
#' \item{for Bayes}{a data.frame with columns: item_id, item_score, mean_beta, SD_beta, -bayes_hpd_b_left-, -bayes_hpd_b_right-}
#' }
#' 
#' 
#' 
coef.prms = function(object, bayes_hpd_b = 0.95, ...)
{

  
  if(bayes_hpd_b <= 0 ||  bayes_hpd_b >= 1)
    stop('args$bayes_hpd_b must be between 0 and 1')
  
  x = object
  
  hpd=function(x, conf=bayes_hpd_b){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(c(l=x[ nnn ],r=x[ n-nn+nnn ]))
  }
  
  if (x$inputs$method=="CML")
  {
    atab=data.frame(item_id=x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    item_score=as.integer(x$inputs$ssIS$item_score[-x$inputs$ssI$first]),
                    beta=x$est$beta,
                    SE_beta=sqrt(diag(x$est$acov.beta)),stringsAsFactors=FALSE)
  }
  
  if (x$inputs$method=="Bayes"){
    
    
    hh = t(apply(x$est$beta,2,hpd))
    atab=data.frame(item_id = x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    a = x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                    mb = colMeans(x$est$beta),
                    sdb = apply(x$est$beta, 2, sd),
                    hpdl = hh[,1], hpdr=hh[,2],stringsAsFactors=FALSE)
    colnames(atab)=c("item_id" ,"item_score", "mean_beta", "SD_beta", 
                     sprintf("%i_hpd_b_left", round(100 * bayes_hpd_b)),
                     sprintf("%i_hpd_b_right", round(100 * bayes_hpd_b)))
  }
  atab$item_id = as.character(atab$item_id)
  rownames(atab) = NULL
  return(atab)
}


#' Functions of theta
#' 
#' returns information function, expected score function or score simulation function 
#' for a single item, an arbitrary group of items or all items
#' 
#' @param parms object produced by fit_enorm or a data.frame with columns item_id, item_score and, 
#' depending on parametrization, a column named beta/delta, eta or b
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
#' \item{var_score}{conditional score variance}
#' }
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
#' plot(density(pv$PV1), col='green', axes=FALSE,xlab=NA, ylab=NA,main=NA)
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
variance_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='variance')
}


theta_function = function(parms, items=NULL, booklet=NULL, which.draw=NULL, 
                          what=c('information','expected','sim', 'variance'))
{
  what = match.arg(what)
  
  # data preparation
  # create fl(item_id,first,last), a, b
  
  if(inherits(parms,'data.frame'))
  {
    out = transform.df.parms(parms,'b',TRUE)
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
    # to do: check if mst needs adjustment for this to work
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
  }  
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
  } else if (what=="variance")
  {
    max_score = sum(a[fl$last])
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta))) 
        stop('theta may not contain nan/NA values') 
      
      distr = pscore(theta,  b=b, a=a, first=fl$first, last=fl$last)
      res = apply(distr, 2, function(x) sum((0:max_score)^2*x)-sum((0:max_score)*x)^2)
      res
    }
    class(out) = append('var_func',class(out))
  }
  
  out
}

print.inf_func = function(x,...) cat('Information function I(theta)\n')
print.exp_func = function(x,...) cat('Conditional expected score function E(X|theta)\n')
print.sim_func = function(x,...)  cat('(x_i1, ..., x_in) ~ ENORM([theta])\n\tfunction (theta)\n')
print.var_func = function(x,...) cat('Conditional score variance var(X|theta)\n')

