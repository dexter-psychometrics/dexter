
#' Latent correlations
#'
#' Estimates correlations between latent traits. Use an item_property to distinguish the different scales. 
#' This function uses plausible values so results may differ slightly between calls. 
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
  # why not check which is the case?
  
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
               MoreArgs = list(method="EAP", prior="Jeffreys"))
  
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
  # parms = cml
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
  
  max_cores = get_ncores(desired = 256, maintain_free = 1L)
  for (i in 1:nIter)
  {
    for (d in 1:nd)
    {
      cons = condMoments(prior$mu, prior$Sigma, d, x.value=pv[,-d]) 
      
      PV_sve(models[[d]]$b, models[[d]]$a, models[[d]]$design$first, models[[d]]$design$last, 					
             models[[d]]$bcni,
             respData[[d]]$x$booklet_id, respData[[d]]$x$booklet_score, cons$mu, sqrt(cons$sigma),
             max_cores,
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
  
  res = lapply(res, function(x){colnames(x) = rownames(x) = names(models); x})
  res$n_persons = np
  res
}


# mean and variance of Y|x if x,y is multivariate normal
#
# with mean mu and variance-covariance matrix sigma
# @param m vector of means
# @param sigma covariance matrix
# @param y.ind indices dependent variable(s)
# @param x.ind indices conditioning variables. If null its just all others
# @param x.value value of conditioning variables
# 
condMoments = function(mu, sigma, y.ind, x.ind=NULL, x.value )
{
  if (is.null(x.ind)) x.ind = setdiff(1:length(mu), y.ind)
  B = sigma[y.ind, y.ind]
  C = sigma[y.ind, x.ind, drop = FALSE]
  D = sigma[x.ind, x.ind]
  CDinv = C %*% solve(D)
  if (is.vector(x.value))
  {
    cMu = c(mu[y.ind] + CDinv %*% (x.value - mu[x.ind]))
  }else
  {
    nP = nrow(x.value)
    cMu = rep(0,nP)
    for (i in 1:nP) cMu[i] = mu[y.ind] + CDinv %*% (x.value[i,] - mu[x.ind])
  }
  cVar = B - CDinv %*% t(C)
  return(list(mu=cMu, sigma=cVar))
}

update_MVNprior = function(pvs,Sigma)
{
  m_pv = colMeans(pvs)
  nP = nrow(pvs)
  mu = rmvnorm(1, mean=m_pv,sigma=Sigma/nP)
  
  S=(t(pvs)-m_pv)%*%t(t(pvs)-m_pv)
  Sigma = solve(rWishart(1,nP-1,solve(S))[,,1]) 
  return(list(mu = mu, Sigma = Sigma))
}


# borrowed the following from mvtnorm source for the time being
# so we don't have to import a whole package
rmvnorm = function(n,mean,sigma)
{
  ev = eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  retval = matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% R
  retval = sweep(retval, 2, mean, "+")
  colnames(retval) = names(mean)
  retval
}

