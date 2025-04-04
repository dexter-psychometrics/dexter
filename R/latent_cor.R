
#' Latent correlations
#'
#' Estimates correlations between latent traits using plausible values as described in Marsman, et al. (2022). 
#' An item_property is used to distinguish the different scales. 
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
#' @return List containing a estimated correlation matrix, the corresponding standard deviations, 
#' and the lower and upper limits of the highest posterior density interval and the complete mcmc sample
#' @details
#' This function uses plausible values so results may differ slightly between calls. 
#' 
#' @references 
#' Marsman, M., Bechger, T. M., & Maris, G. K. (2022). Composition algorithms for conditional distributions. 
#' In Essays on Contemporary Psychometrics (pp. 219-250). Cham: Springer International Publishing.
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
  
  from = 50L
  by = 2L
  which.keep = seq(from,(from-by)*(from>by)+by*nDraws,by=by)
  nIter = max(which.keep)
  
  if(is.matrix(dataSrc))
  {
    stop("a matrix datasrc is not yet implemented for this function")
  }
  # merge_within_persons might also be user set. But TRUE will be the most useful in nearly all cases
  
  respData = get_resp_data(dataSrc, qtpredicate, env=env, extra_columns=item_property,
    merge_within_persons=TRUE)
  
  respData$x[[item_property]] = ffactor(as.character(respData$x[[item_property]]))
  lvl = levels(respData$x[[item_property]])
  nd = length(lvl)
  
  pb$set_nsteps(nIter + 4*nd)
  
  # this assumes merge over booklets
  persons = respData$x |>
    distinct(.data$person_id, .data[[item_property]]) |>
    count(.data$person_id)
  
  persons = persons |>
    filter(.data$n==nd) |>
    mutate(new_person_id = dense_rank(.data$person_id)) |>
    select('person_id','new_person_id')
  
  np = nrow(persons)
  
  if(np==0) stop_("There are no persons that have scores for every subscale.") 
  
  respData = lapply(split(respData$x, respData$x[[item_property]]), get_resp_data)
  
  models = lapply(respData, function(x){try(fit_enorm(x), silent=TRUE)})
  names(models) = names(respData)
  
  if(any(sapply(models,inherits,what='try-error')))
  {
    message('\nThe model could not be estimated for one or more item properties, reasons:')
    models[!sapply(models,inherits,what='try-error')] = 'OK'
    print(data.frame(item_property=names(models), 
      result=gsub('^.+try\\(\\{ *:','',trimws(as.character(models)))))
    stop('Some models could not be estimated')
  }
  
  for(i in seq_along(respData))
  {
    respData[[i]] = get_resp_data(respData[[i]], summarised=TRUE)
    
    respData[[i]]$x = respData[[i]]$x |> 
      inner_join(persons,by='person_id') |>
      select(-'person_id') |>
      rename(person_id='new_person_id')
  }
  
  pb$tick(2*nd)
  
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
  
  acor = cor(pv)
  # attenuation correction
  for (i in 1:(nd-1))
  {
    for (j in ((i+1):nd)){
      acor[i,j] = acor[i,j]/sqrt(reliab[i]*reliab[j])
      acor[j,i] = acor[i,j]
    }
  }
  # ignore attenuation if result is invalid
  if(any(abs(acor)>1))
    acor = cor(pv)
  
  # make everything simple and zero indexed
  
  models = lapply(models, simplify_parms)
  respData = lapply(respData, get_resp_data, summarised=TRUE)
  
  for(d in 1:nd)
  {
    respData[[d]]$x$booklet_id = as.integer(respData[[d]]$x$booklet_id) - 1L
    models[[d]]$bcni = c(0L,cumsum(table(as.integer(models[[d]]$design$booklet_id))))
  }
  
  prior = list(mu=mean_pv, Sigma = diag(1.7*sd_pv) %*% acor %*% diag(1.7*sd_pv))
  
  mcmc = array(0, dim = c(length(which.keep), nd, nd))
  tel = 1
  
  max_cores = get_ncores(desired = 128L, maintain_free = 1L)
  for (i in 1:nIter)
  {
    for (d in 1:nd)
    {
      cons = condMoments(mu = prior$mu, sigma = prior$Sigma, index=d, value=pv) 
      
      PV_sve(models[[d]]$b, models[[d]]$a, models[[d]]$design$first0, models[[d]]$design$last0, 					
        models[[d]]$bcni,
        respData[[d]]$x$booklet_id, respData[[d]]$x$booklet_score, cons$mu, sqrt(cons$sigma),
        max_cores,
        pv,d-1L, 10L)
    }
    
    prior = update_MVNprior(pv,prior$Sigma)
    
    if (i %in% which.keep){
      mcmc[tel,,] = cov2cor(prior$Sigma)
      tel = tel+1L
    }
    pb$tick()
  }
  
  hp = apply(mcmc,2:3,hpdens)
  
  res = list(cor = apply(mcmc,2:3,mean), sd = apply(mcmc,2:3,sd),
    hpd_l = hp[1,,], hpd_h = hp[2,,])
  
  res = lapply(res, function(x){colnames(x) = rownames(x) = names(models); x})
  res$n_persons = np
  res$mcmc = mcmc
  res
}

# robust for not entirely positive definite matrix V
update_MVNprior = function(pvs,Sigma)
{
  m_pv = colMeans(pvs)
  nP = nrow(pvs)
  mu = rmvnorm(1, mu=m_pv, sigma=Sigma/nP)
  
  pvs = sapply(1:ncol(pvs),\(i) pvs[,i] - m_pv[i])
  S = t(pvs) %*% pvs

  V = solve(S)
  
  Sigma = try({solve(rWishart(1,nP-1,V)[,,1])}, silent=TRUE)
  if(inherits(Sigma,'try-error'))
  {
    V = nearPD(V)
    Sigma = solve(rWishart(1,nP-1,V)[,,1])
  }   
  return(list(mu = mu, Sigma = Sigma))
}

# mean and variance of Y given values for x
# near zero/negative variances can occur,  meaningfull error message?

condMoments = function(mu, sigma, index, value )
{
  C = sigma[index,-index,drop=FALSE]
  D = sigma[-index,-index]
  CDinv = C %*% solve(D)
  
  mu_y = apply(value[,-index,drop=FALSE],1, `-`, mu[-index])
  
  if(is.null(dim(mu_y)))
  {
    mu_y = (CDinv %*% mu_y) + mu[index]
  } else
  {
    mu_y = apply(mu_y,2,\(x) CDinv %*% x) + mu[index]
  }

  return(list(mu=mu_y, sigma = sigma[index,index] - CDinv %*% t(C)))
}


# somewhat simplified, from the Matrix package
# Jens Oehlschlägel and Matrix package authors
# used in edge case in pudate mvn prior
nearPD = function (x,  eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxit = 50L) 
{
  if (!isSymmetric(x)) {
    x = (x + t(x))/2 
  }
  
  n = ncol(x)
  
  D_S = x
  D_S[] = 0
  
  X = x
  iter = 0L
  converged = FALSE
  conv = Inf
  while (iter < maxit && !converged) {
    Y = X
    
    R = Y - D_S
    e = eigen(R)
    Q = e$vectors
    d = e$values
    p = d > eig.tol * d[1]
    if (!any(p)) 
      stop("Matrix seems negative semi-definite")
    Q = Q[, p, drop = FALSE]
    X = tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
    
    D_S = X - R
    
    conv = norm(Y - X, "I")/norm(Y, "I")
    iter = iter + 1L
    
    converged = (conv <= conv.tol)
  }
  
  e = eigen(X, symmetric = TRUE)
  d = e$values
  Eps = posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] = Eps
    
    Q = e$vectors
    o.diag = diag(X)
    X = Q %*% (d * t(Q))
    D = sqrt(pmax(Eps, o.diag)/diag(X))
    X[] = D * X * rep(D, each = n)
  }
  
  X
}


rmvnorm = function(n, mu, sigma)
{
  ev = eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  res = matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% R
  res = sweep(res, 2, mu, "+")
  colnames(res) = names(mu)
  res
}
