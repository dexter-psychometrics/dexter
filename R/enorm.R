
##########################################
#' Fit the extended nominal response model
#'
#' Fits the extended nominal response model by Bayesian sampling or CML
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param fixed_params Optionally, a prms object from a previous analysis or 
#' a data.frame with columns: item_id, item_score (omitting 0 score category) and beta
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nIterations Number of Gibbs samples when method is Bayes, max. number of iterations 
#' when method is CML
#' @return An object of type \code{prms}. The prms object can be cast to a data.frame of item parameters 
#' using the \code{as.data.frame} built-in function or used directly as input for other Dexter functions.
#' 
#' @seealso functions that accept a prms object as input: \code{\link{ability}}, \code{\link{plausible_values}}
#'
fit_enorm = function(dataSrc, predicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), nIterations=500)
{
  method = match.arg(method)
  qtpredicate = eval(substitute(quote(predicate)))
  fit_enorm_(dataSrc, qtpredicate = qtpredicate, fixed_params = fixed_params, method=method, nIterations=nIterations, env=caller_env())
}


fit_enorm_ = function(dataSrc, qtpredicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), nIterations=500, env=NULL) 
{
  method <- match.arg(method)
  if(is.null(env)) env = caller_env()
  
  r = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env)
  
  x = r$x
  design = r$design
  if(nrow(x) == 0) stop('no data to analyse')
  if(nrow(design) == 1) stop('There are responses to only one item in your selection, this cannot be calibrated.') 
  
  if(length(unique(design$booklet_id)) > 1)
  {
    im = as.matrix(table(design$item_id, design$booklet_id))
    wm = crossprod(im, im)
    diag(wm) = 0
    if(!design_is_connected(list(im=im, wm=wm))) stop('Your design is not connected')  
  }
  
  itm_max = x %>% 
    group_by(.data$item_id) %>% 
    summarise(maxScore=max(.data$item_score)) %>%
	  ungroup()
  
  if(any(itm_max$maxScore==0)) stop('One or more items has a maximum score of 0')
  # for now, just error. May possibly become warning later

  ssBIS = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$item_score) %>% 
    summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>% 
    ungroup()
  
  plt = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$sumScore) %>% 
    summarise(meanScore=mean(.data$item_score), N=n()) %>% 
    ungroup()
  
  maxScores = itm_max %>%
    inner_join(design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    summarise(maxTotScore = sum(.data$maxScore))
  
  allScores = maxScores %>% group_by(.data$booklet_id) %>%
    do(tibble(sumScore=0:.$maxTotScore))
  
  stb = plt %>%
    select(.data$booklet_id, .data$sumScore, .data$N) %>%
    distinct() %>%
    right_join(allScores, by=c('booklet_id','sumScore')) %>%
    do({
      .$N[is.na(.$N)]=0
      .$booklet_id[is.na(.$booklet_id)] = .$booklet_id[!is.na(.$booklet_id)][1]
      as.data.frame(.)
    }) %>%
    select(.data$booklet_id, .data$sumScore, .data$N) %>%
    arrange(.data$booklet_id, .data$sumScore)
  
  ssIS = ssBIS %>% 
    group_by(.data$item_id, .data$item_score) %>%
    summarise(sufI = sum(.data$sufI)) %>%
    ungroup() %>%
    arrange(.data$item_id, .data$item_score)
  
  ssI  = ssIS %>% 
    group_by(.data$item_id) %>%
    summarise(nCat = n()) %>% 
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1,last = cumsum(.data$nCat)) %>%
    ungroup() %>%
    arrange(.data$item_id)
  
  if(any(ssI$nCat == 1)) stop('One or more items are without score variation')
  
  m = x  %>% 
    group_by(.data$booklet_id)  %>% 
    summarise(m = n_distinct(.data$person_id))
  
  a = ssIS$item_score
  
  
  bkl = lapply(m$booklet_id, function(x) 
  {
    itInBk = design$item_id[design$booklet_id==x] 
    items = ssI$item_id[ssI$item_id %in% itInBk]
    m = m$m[m$booklet_id==x]
    first = ssI$first[match(items, ssI$item_id)]
    last =  ssI$last[match(items, ssI$item_id)]
    list(booklet=x, items=items, m=m, first=first, last=last, 
         score=stb$sumScore[stb$booklet_id==x], 
         scoretab=stb$N[stb$booklet_id==x], 
         lambda=rep(1,sum(a[last])+1))
  })
  
  names(bkl) = m$booklet_id
  
  itemList = lapply(ssI$item_id, function(x) design$booklet_id[design$item_id==x])
  
  itemListInt = lapply(itemList, function(x)match(x,m$booklet_id))
  
  it_sc_lab = paste0(ssIS$item_id[-ssI$first], "_",ssIS$item_score[-ssI$first])
  
  ## deal with fixed parameters
  if(is.null(fixed_params))
  { 
    fixed_b = NULL
  } else 
  {  
    if(inherits(fixed_params,'prms'))
    {
      fixed_b = (
        fixed_params$inputs$ssIS %>%
          bind_cols(fixed_params$est$b) %>%
          right_join(ssIS, by=c('item_id','item_score')) %>%
          arrange(.data$item_id,.data$item_score)
      )$b
    } else
    {
      # transform the oplmlike fixed params to the parametrization dexter uses internally
      
      #some cleaning and checking
      fixed_params = fixed_params %>% mutate_if(is.factor, as.character) 
      
      if(length(setdiff(c('item_id','item_score','beta'),colnames(fixed_params))) > 0)
      {
        stop(paste('fixed_params does not contain required column(s):',
                   paste0(setdiff(c('item_id','item_score','beta'),colnames(fixed_params)),collapse=',')))
      }
      # check for missing categories in fixed_params
      # missing categories in data are presumably not a problem
      ssIS %>% 
        semi_join(fixed_params, by='item_id') %>%
        left_join(fixed_params, by=c('item_id','item_score')) %>%
        filter(is.na(.data$beta) & .data$item_score !=0) %>%
        do({
          if(nrow(.) > 0)
          {
            print('The following score categories are missing:')
            print(as.data.frame(
              .[,c('item_id','item_score')] %>% 
                  arrange(.data$item_id, .data$item_score)))
            stop('missing score categories in fixed_params, see output')
          }
          data.frame()
        })
      
      # we have to mimic the positional dexter structure to be able to call toDexter
      # omit itemscores not found in data
      # an alternative would be to keep them, but that would possibly mess up the prms object
      # I do wondr if it is better to throw them out before or after the toDexter call
      fixed_params = fixed_params %>%
        semi_join(ssIS, by=c('item_id','item_score')) %>%
        filter(.data$item_score != 0) %>%
        arrange(.data$item_id, .data$item_score) 
      
      fixed_ssI  = fixed_params %>% 
        group_by(.data$item_id) %>%
        summarise(nCat = n()) %>% 
        mutate(first = cumsum(.data$nCat) - .data$nCat + 1,last = cumsum(.data$nCat)) %>%
        ungroup() %>%
        arrange(.data$item_id)
      
      # now we will unfortunately loose the item id's 
      # which we'll need to join the fixed and unfixed again to get them in the right order
      # we also have to include the 0 category again
      dx_b = toDexter(fixed_params$beta, fixed_params$item_score, fixed_ssI$first, fixed_ssI$last, 
                      re_normalize=FALSE)
      dx_b = tibble(item_score = dx_b$inputs$ssIS$item_score, 
                    b = dx_b$est$b,
                    item_id = (fixed_params %>%
                                 group_by(.data$item_id) %>% 
                                 do({tibble(item_id=c(.$item_id,as.vector(.$item_id)[1]))})
                                )$item_id )

      # make a tibble of the fixed items including the 0 categories and 
      # cbind with the dexter parametrization
      # then join with all the items again
      fixed_b = dx_b %>%
        full_join(ssIS, by=c('item_id','item_score')) %>%
        arrange(.data$item_id, .data$item_score) %>%
        pull(.data$b) 
    }
    if(!any(is.na(fixed_b))) stop('nothing to calibrate, all parameters are fixed')
  }
  
  if (method=="CML"){
    result = try(calibrate_CML(booklet=bkl, sufI=ssIS$sufI, a=ssIS$item_score, 
                               first=ssI$first, last=ssI$last, nIter=nIterations,
                               fixed_b=fixed_b))
    names(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score)
    row.names(result$beta.cml)=it_sc_lab
    rownames(result$acov.cml)=it_sc_lab
    colnames(result$acov.cml)=it_sc_lab
  } else {

    b = exp(runif(nrow(ssIS), -1, 1))
    result = try(calibrate_Bayes(itemList=itemListInt, booklet=bkl, sufI=ssIS$sufI, b=b, a=a, 
                                 first=ssI$first, last=ssI$last, nIter=nIterations,
                                 fixed_b=fixed_b))
    colnames(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score) 
    colnames(result$beta.cml)=it_sc_lab
  }
  if (inherits(result, "try-error")) result=NULL
  
  outpt = list(est=result, 
               inputs=list(bkList=bkl, ssIS=ssIS, ssI=ssI, 
                           stb=stb, method=method), 
               xpr=as.character(qtpredicate))
  class(outpt) = append('prms', class(outpt)) 
  outpt
}




##########################################
#' A print method for ENORM parameters
#'
#' @param x An object produced by function \code{fit_enorm}
#' @param ... Any other parameters to the print method
#' @method print prms
#'
#'
print.prms <- function(x, ...){
  
  hpd=function(x, conf=0.95){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(paste0("(",as.character(round(x[ nnn ],digits=3))," , "
                  ,as.character(round(x[ n-nn+nnn ],digits=3)),")") )
  }
  
  if (x$inputs$method=="CML")
  {
    atab=as.data.frame(cbind(x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                             x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                             round(x$est$beta.cml,digits=3),
                             round(sqrt(diag(x$est$acov.cml)),digits=3)))
    colnames(atab)=c("item_id" ,"a", "Beta", "SE(Beta)")
  }
  
  if (x$inputs$method=="Bayes"){
    hh=apply(x$est$beta.cml,2,hpd)
    atab=as.data.frame(cbind(x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                             x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                             round(colMeans(x$est$beta.cml),digits=3),
                             round(apply(x$est$beta.cml, 2, sd),digits=3),
                             hh))
    colnames(atab)=c("item_id" ,"a", "mean(Beta)", "SD(Beta)", "95%hpd(B)")
  }
  row.names(atab)=NULL
  do.call(print,append(list(x=atab),list(...)))
}

##########################################
#' Coerce enorm parameters object to a data.frame of parameters
#'
#' @param x An object produced by function \code{fit_enorm}
#' @param row.names	NULL or a character vector giving the row names for the data frame. 
#' Missing values are not allowed.
#' @param optional Will be set to FALSE, so ignore.
#' @param ... any other parameters to the as.data.frame method are ignored silently
#' @method as.data.frame prms
#'
#'
as.data.frame.prms <- function(x, row.names=NULL, optional=FALSE, ...){

  hpd=function(x, conf=0.95){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(data.frame(l=round(x[ nnn ],digits=3),r=round(x[ n-nn+nnn ],digits=3)))
  }
  
  if(optional) optional=FALSE
  if (x$inputs$method=="CML")
  {
    atab=data.frame(item_id=x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    a=as.integer(x$inputs$ssIS$item_score[-x$inputs$ssI$first]),
                    B=x$est$beta.cml,
                    se=sqrt(diag(x$est$acov.cml)),stringsAsFactors=FALSE)
    colnames(atab)=c("item_id" ,"item_score", "beta", "SE_b")
  }
  
  if (x$inputs$method=="Bayes"){
    hh = map_df(apply(x$est$beta.cml,2,hpd),c)
    atab=data.frame(item_id = x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    a = x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                    mb = colMeans(x$est$beta.cml),
                    sdb = apply(x$est$beta.cml, 2, sd),
                    hpdl = hh[,1], hpdr=hh[,2],stringsAsFactors=FALSE)
    colnames(atab)=c("item_id" ,"item_score", "mean_beta", "SD_b", "95_hpd_b_left","95_hpd_b_right")
    
  }
  rownames(atab) = row.names
  return(atab)
}




##########################################
#' Estimate abilities
#'
#' Computes estimates of ability and optionally attaches them to person data by sum score
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param method If ability will be estimated with maximum likelihood (MLE). Otherwise, 
#' we produce a Bayesian expected a posteriori (EAP) estimate.
#' @param use_draw When the ENORM was fitted with a Gibbs sampler (this is 
#' recognised automatically), the number of the random draw (iteration) to use. 
#' If NULL, all draws will be averaged. If outside range,
#' the last iteration will be used.
#' @param person_level If TRUE, return results per person, otherwise just
#' the score transformation table.   
#' @return if \code{person_level=TRUE} a data.frame with columns:
#' (booklet_id, person_id, sumScore, theta), if \code{person_level=FALSE}
#' a data.frame with columns: (booklet_id, sumScore, theta) 
#'
#' @details MLE estimates of ability will produce an NA for
#' the minimum (=0) or the maximum score on a booklet. If this is undesirable, 
#' we advise to use the EAP estimate.
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' f = fit_enorm(db)
#' aa = ability(db,f,method="MLE",person_level = FALSE)
#' bb = ability(db,f,method="EAP",person_level = FALSE)
#' plot(bb$sumScore, bb$theta, xlab="test-score", ylab="ability est.", pch=19, cex=0.7)
#' points(aa$sumScore, aa$theta, col="red", pch=19, cex=0.7)
#' legend("topleft", legend = c("EAP", "MLE"), bty = "n",
#'     lwd = 1, cex = 0.7, col = c("black", "red"), lty=c(0,0), pch = c(19,19))
#' 
#' close_project(db)
#' }
#' 
#' 
ability <- function(dataSrc, parms, predicate=NULL, method=c("MLE","EAP"), use_draw=NULL, person_level=TRUE){
  
  method <- match.arg(method)
  qtpredicate=eval(substitute(quote(predicate)))
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=caller_env())

  if(nrow(respData$x)==0) stop('no data to analyse')
  
  # we make the design and join with the indexes
  design = respData$design %>%
    left_join(parms$inputs$ssI, by='item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  if(any(is.na(design$first))) stop('some of your items are without parameters')
  
  # check if all scores are known
  
  if(nrow(anti_join(respData$x, parms$inputs$ssIS, by=c('item_id','item_score'))) > 0)
  {
    stop('Some item_scores in your data are not present in your parameters.')
  } 
  # now we can summarise
  respData = get_resp_data(respData, summarised=TRUE)
  
  if(parms$input$method=="CML"){
    b = parms$est$b
    a = parms$inputs$ssIS$item_score
  } else {
    a = parms$est$a
    if(is.null(use_draw)) {
      b = colMeans(parms$est$b)  
    } else {
      if (use_draw %in% 1:nrow(parms$est$b)) {
        b = parms$est$b[use_draw,]   
      } else {
        b = parms$est$b[nrow(parms$est$b),]
      }
    }   
  } 
  
  mx_scores = parms$inputs$ssIS %>%  
    group_by(.data$item_id) %>% 
    summarise(item_max = max(.data$item_score))
  
  ou = design %>% 
    inner_join(mx_scores, by='item_id') %>% 
    group_by(.data$booklet_id) %>% 
    summarise(maxScore = sum(.data$item_max))
  
  ou = merge(ou, data.frame(sumScore=0:max(ou$maxScore))) %>%
    filter(.data$sumScore <= .data$maxScore) %>%
    select(.data$booklet_id, .data$sumScore) %>%
    arrange(.data$booklet_id, .data$sumScore)
  
  
  if (method=="MLE")
  {
    ou$theta = unlist(map(sort(unique(ou$booklet_id)), 
                          ~theta_MLE(b,a,design[design$booklet_id==.x,]$first, design[design$booklet_id==.x,]$last)))
  }else
  {
    ou$theta = unlist(map(sort(unique(ou$booklet_id)),
                          ~theta_EAP(b,a,design[design$booklet_id==.x,]$first, design[design$booklet_id==.x,]$last,
                                     ou$sumScore[ou$booklet_id==.x])))    
  }
  
  if (!person_level)
  {
    return(ou) 
  } else{
    return(respData$x %>% 
             inner_join(ou, by = c("booklet_id", "sumScore")) %>% 
             select(.data$booklet_id, .data$person_id, .data$sumScore, .data$theta))
  }
}

