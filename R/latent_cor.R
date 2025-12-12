

#' Latent correlations
#'
#' Estimates correlations between latent traits using plausible values as described in Marsman, et al. (2022). 
#' An item_property is used to distinguish the different scales. 
#'
#' @param dataSrc A connection to a dexter database or a data.frame with columns: person_id, item_id, item_score and 
#' the item_property
#' @param item_property The name of the item property used to define the domains. If \code{dataSrc} is a dexter db then the
#' item_property must match a known item property. If dataSrc is a data.frame, item_property must be equal to
#'  one of its column names.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param nDraws Number of draws for plausible values
#' @param hpd deprecated, use the `coef` method to set the highest posterior density interval.
#' @param use complete.obs uses only persons with answers on all domains. Pairwise.complete.obs uses all cases 
#' for which there are responses in at least two domains.
#' @return `latent_cor` object, which is a list containing an estimated (mean) correlation matrix, the corresponding standard deviations, 
#' and the complete mcmc sample. Use the coef method to extract highest posterior density intervals around the estimated correlation matrix.
#' @details
#' To compute latent correlations, a model is estimated for each subscale. If a design for any subscale is not connected this will result in an error.  
#' 
#' Latent correlations are generated using a Bayesian approach. `Use` is `"pairwise.complete.obs"` works slightly different from 
#' `cor` since complete matrices are imputed. Therefore individual correlation matrices in the mcmc sample are positive semi-definite even with the pairwise option 
#' (assuming the matrix is not degenerate). However the mean of the mcmc sample need never be positive semidefinite. 
#'  
#' @seealso \code{\link{coef.latent_cor}}
#' 
#' @references 
#' Marsman, M., Bechger, T. M., & Maris, G. K. (2022). Composition algorithms for conditional distributions. 
#' In Essays on Contemporary Psychometrics (pp. 219-250). Cham: Springer International Publishing.
#' 
latent_cor = function(dataSrc, item_property, predicate=NULL, nDraws=500, hpd=0.95, use=c("complete.obs","pairwise.complete.obs"))
{
  
  check_dataSrc(dataSrc, matrix_ok=FALSE)
  check_num(nDraws, 'integer', .length=1, .min=1)
  check_num(hpd,  .length=1, .min=1/nDraws,.max=1)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  use = match.arg(use)
  if(!missing(hpd))
    cl_msg('Argument `hpd` is deprecated, use %s to set a confidence envelope for the highest posterior density interval.', 
      'coef.latent_cor', mod = 'message')


  pb = get_prog_bar(nDraws, retrieve_data = is_db(dataSrc), lock=TRUE)
  on.exit({pb$close()})
  
  from = 50L
  by = 2L
  which.keep = seq(from,(from-by)*(from>by)+by*nDraws,by=by)
  nIter = max(which.keep)
  
  # merge_within_persons might also be user set. But TRUE will be the most useful in nearly all cases
  respData = get_resp_data(dataSrc, qtpredicate, env=env, extra_columns=item_property,
    merge_within_persons=TRUE)
  
  respData$x[[item_property]] = ffactor(as.character(respData$x[[item_property]]))
  ndomains = nlevels(respData$x[[item_property]])
  
  if(ndomains == 1)
    stop_('`item_property` should define at least two domains')
  
  pb$set_nsteps(nIter + 4*ndomains)
  
  # this assumes merge over booklets is TRUE
  persons = respData$x |>
    distinct(.data$person_id, .data[[item_property]]) |>
    count(.data$person_id) |>
    filter(.data$n >= if(use=="complete.obs") ndomains else 2L) |>
    mutate(new_person_id = dense_rank(.data$person_id)) |>
    select('person_id','new_person_id')

  np = nrow(persons)
  if(np == 0)
      stop_("There are no persons that have scores for every subscale.")


  respData = lapply(split(respData$x, respData$x[[item_property]]), get_resp_data)
  
  models = lapply(respData, function(x){try(fit_enorm(x), silent=TRUE)})
  names(models) = names(respData)
  
  if(any(sapply(models,inherits,what='try-error')))
  {
    message('\nThe model could not be estimated for one or more item properties, reasons:')
    models[!sapply(models,inherits,what='try-error')] = 'OK'
    print(data.frame(item_property=names(models), 
      result=gsub('^.+try\\(\\{ *:','',trimws(as.character(models)))))
    stop_('Some models could not be estimated')
  }
  
  
  for(i in seq_along(respData))
  {
    respData[[i]] = get_resp_data(respData[[i]], summarised=TRUE)
    
    respData[[i]]$x = respData[[i]]$x |> 
      inner_join(persons,by='person_id') |>
      select(-'person_id') |>
      rename(person_id='new_person_id')
  }
  
  pb$tick(2*ndomains)
  
  pv = matrix(NA_real_,np,ndomains)
  
  reliab = rep(0,ndomains)
  sd_pv = rep(0,ndomains)
  mean_pv = rep(0,ndomains)
  
  for (i in 1:ndomains)
  {
    pvs = plausible_values(respData[[i]],models[[i]],nPV = 2)
    reliab[i] = cor(pvs$PV1,pvs$PV2)
    sd_pv[i] = sd(pvs$PV1)
    mean_pv[i] = mean(pvs$PV1)
    pv[pvs$person_id,i] = pvs$PV1
  }
  pb$tick(2*ndomains)
  
  # If NA entries occur they will be filled with a guesstimate
  acor = raw_cor = cor_fill_na(cor(pv, use = "pairwise.complete.obs"))
  
  # attenuation correction
  for (i in 1:(ndomains-1))
  {
    for (j in ((i+1):ndomains)){
      acor[i,j] = acor[i,j]/sqrt(reliab[i]*reliab[j])
      acor[j,i] = acor[i,j]
    }
  }
  # reverse attenuation if result is invalid
  if(any(abs(acor)>1))
    acor = raw_cor
  
  # make everything simple and zero indexed
  
  models = lapply(models, simplify_parms)

  for(dom in 1:ndomains)
  {
    respData[[dom]]$x$booklet_id = as.integer(respData[[dom]]$x$booklet_id) - 1L
    models[[dom]]$bcni = c(0L,cumsum(table(as.integer(models[[dom]]$design$booklet_id))))
  }
  
  prior = list(mu=mean_pv, Sigma = diag(1.7*sd_pv) %*% acor %*% diag(1.7*sd_pv))
  
  mcmc = array(0, dim = c(length(which.keep), ndomains, ndomains))
  tel = 1
  missing_data = if(use=='complete.obs') NULL else apply(pv,2,is.na) + 0L

  #impute missing
  pv = mice(pv,acor) 
  max_cores = get_ncores(desired = 128L, maintain_free = 1L)


  for (i in 1:nIter)
  {
    for (dom in 1:ndomains)
    {
      cons = condMoments(mu = prior$mu, sigma = prior$Sigma, index=dom, value=pv) 
      
      PV_sve(models[[dom]]$b, models[[dom]]$a, models[[dom]]$design$first0, models[[dom]]$design$last0,
       models[[dom]]$bcni,
       respData[[dom]]$x$booklet_id, respData[[dom]]$x$booklet_score, cons$mu, sqrt(cons$sigma),
       max_cores=max_cores,
       pv, missing_data = missing_data[,dom], pv_col_indx = dom-1L, niter=10L)
    }

    prior = update_MVNprior(pv,prior$Sigma)
    
    if (i %in% which.keep){
      mcmc[tel,,] = cov2cor(prior$Sigma)
      tel = tel+1L
    }
    pb$tick()
  }
  
  hp = apply(mcmc,2:3,hpdens,conf=hpd)
  
  res = list(cor = apply(mcmc,2:3,mean), sd = apply(mcmc,2:3,sd),
    hpd_l = hp[1,,], hpd_h = hp[2,,])
  
  res = lapply(res, function(x){colnames(x) = rownames(x) = names(models); x})
  res$n_persons = np
  res$mcmc = mcmc
  res$use = use
  class(res) = append('latent_cor',class(res))
  res
}

print.latent_cor = function(x,...)
{
  d = nrow(x$cor)
  nms = dimnames(x$cor)
  hp = apply(x$mcmc,2:3,hpdens,conf=0.95)
  
  cat(sprintf('n=%i, %s\n\ncorrelation:\n\n',x$n_persons,x$use))
  print(matrix(sprintf('%.2f',x$cor),nrow=d,dimnames=nms),quote=FALSE)
  
  cat('\nhighest posterior density (0.95):\n\n')
  h = matrix(sprintf('%.2f-%.2f',hp[1,,],hp[2,,]),nrow=d,dimnames=nms)
  diag(h)=''
  print(h,quote=FALSE)
}

#' latent correlations
#'
#' @param object resulting value from latent_cor
#' @param hpd width of the confidence envelope of the  Bayesian highest posterior density interval around the correlations, 
#' value must be between 0 and 1.
#' @param ... ignored
#' @returns list with matrices cor, sd, hpd_l (lower boundary), hpd_u (upper boundary), array mcmc with the complete mcmc sample 
#' 
#' Extract correlations with a specified highest posterior density interval
#'  
coef.latent_cor = function(object, hpd=0.95, ...)
{
  hp = apply(object$mcmc,2:3,hpdens,conf=hpd)

  object$hpd_l = hp[1,,]
  object$hpd_h = hp[2,,]
  #object$mcmc = NULL
  class(object) = 'list'
  object
}


# helper functions for when stepwise is used and NA's occur in starting values

# fails if matrix is not connected I think, which is exactly what should happen
cor_fill_na = function(x)
{
  if(!anyNA(x)) return(x)
  
  m = which(is.na(x), arr.ind=TRUE)
  m = m[m[,2] >m[,1],,drop=FALSE]
  
  dup = \(i) sapply(i,\(j) sum(i==j)>1)
  
  # one by one because lbfgsb method often gives weird results
  while(nrow(m) > 1)
  {
    keep = m[!dup(m[,1]),1][1]
    if(is.na(keep) || length(keep)==0) 
    {
      keep = m[!dup(m[,2]),2][1]
      if(is.na(keep) || length(keep)==0) 
        break
      keep = sort(c(keep,setdiff(1:nrow(x),m[,2])))
    } else
    {
      keep = sort(c(keep,setdiff(1:nrow(x),m[,1])))
    }
    
    if(length(keep)<2) break
    x[keep,keep] = cor_fill_na(x[keep,keep])
    m = which(is.na(x), arr.ind=TRUE)
    m = m[m[,2] >m[,1],,drop=FALSE]
  }
  
  start = rep(mean(x[upper.tri(x, diag=FALSE)],na.rm=TRUE),nrow(m))
  
  fn = function(pars)
  {
    x[m] = pars
    x[lower.tri(x)] = t(x)[lower.tri(x)]
    -det(x)
  }
  
  opt = optim(start, fn = fn, method = if(nrow(m)>1) "L-BFGS-B" else "Brent", lower=-.99, upper=.99)
  
  pars = opt$par
  x[m] = pars
  x[lower.tri(x)] = t(x)[lower.tri(x)]
  x
}


# x = rnorm(20000)
# dat = cbind(x,x+rnorm(20000,0,0.2),x+rnorm(20000,0,0.3))
# cor(dat)
# dat[1:10000,2] = NA
# r = cor(dat,use='pair')
# y = mice(dat,r)
mice = function(x,r)
{
  if(!anyNA(x)) return(x)
  
  m = colMeans(x,na.rm=TRUE)
  s = apply(x,2,sd,na.rm=TRUE)
  
  # z scores
  z = sapply(1:ncol(x),\(i) (x[,i]-m[i])/s[i])
  
  for(j in 1:ncol(x))
  {
    i = which(is.na(x[,j]))
    if(length(i)>0)
      x[i,j] = m[j] + s[j] * coalesce(rowMeans(sweep(z[i,-j,drop=FALSE],2,r[j,-j,drop=TRUE],`*`),na.rm=TRUE),0)
  } 

  x  
}

# priors and moments for latent_cor

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


# somewhat simplified, from the Matrix package included with R
# original function by Jens Oehlschlägel and Matrix package authors
# used in edge case in update mvn prior to give some slight relief, but usually when matrices become degenerative
# it will fail later anyhow
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
