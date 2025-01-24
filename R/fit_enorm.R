
#' Fit the extended nominal response model
#'
#' Fits an Extended NOminal Response Model (ENORM) using conditional maximum likelihood (CML)
#' or a Gibbs sampler for Bayesian estimation. 
#' 
#' @details
#' 
#' The eNRM is a slight generalization of the PCM and the OPLM. It
#' reduces to the Rach model for dichotomous items when all itemscores are 0 or 1, is equal to the PCM for polytomous items if all
#' itemscores up to the maximum score occur. It is equal to the oplm if all itemscores have an equal common divisor larger than 1.
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param fixed_params Optionally, an \code{enorm} object from a previous analysis or 
#' a data.frame with parameters, see details.
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nDraws Number of Gibbs samples when estimation method is Bayes. 
#' @param merge_within_persons whether to merge different booklets administered to the same person, enabling linking over persons as well as booklets.
#' @return An object of type \code{enorm}. The enorm object can be cast to a data.frame of item parameters 
#' using function \code{coef} or used directly as input for other Dexter functions.
#' @details
#' To support some flexibility in fixing parameters, fixed_params can be a dexter enorm object or a data.frame.
#' If a data.frame, it should contain the columns item_id, item_score and a difficulty parameter beta
#' 
#' @references 
#' Maris, G., Bechger, T.M. and San-Martin, E. (2015) A Gibbs sampler for the (extended) marginal Rasch model. 
#' Psychometrika. 80(4), 859-879. 
#' 
#' Koops, J. and Bechger, T.M. and Maris, G. (2024); Bayesian inference for multistage and other 
#' incomplete designs. In Research for Practical Issues and Solutions in Computerized Multistage Testing.
#' Routledge, London. 
#' 
#' @seealso functions that accept an \code{enorm} object as input: \code{\link{ability}}, \code{\link{plausible_values}}, 
#' \code{\link{plot.enorm}}, and \code{\link{plausible_scores}}
#'
fit_enorm = function(dataSrc, predicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), 
                     nDraws=1000, merge_within_persons=FALSE)
{
  method = match.arg(method)
  qtpredicate = eval(substitute(quote(predicate)))
  env = rlang::caller_env()
  check_dataSrc(dataSrc)

  fit_enorm_(dataSrc, qtpredicate = qtpredicate, fixed_params = fixed_params,
             method=method, nDraws=nDraws, env=env, merge_within_persons=merge_within_persons)
}


fit_enorm_ = function(dataSrc, qtpredicate = NULL, fixed_params = NULL, method=c("CML", "Bayes","Bayes0"), 
                      nDraws=1000, env=NULL, merge_within_persons=FALSE) 
{
  method = match.arg(method)
  if(is.null(env)) env = caller_env()

  pb = get_prog_bar(retrieve_data = is_db(dataSrc), nsteps = if.else(method=='Bayes',nDraws,NULL))
  on.exit({pb$close()})

  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env, retain_person_id=FALSE,
                           merge_within_persons = merge_within_persons)
  
  pb$tick()
  ss =  get_sufStats_nrm(respData)
  
  ## deal with fixed parameters
  ## maybe use simplify parms?
  if(!is.null(fixed_params))
  {
    if(inherits(fixed_params,'enorm') || inherits(fixed_params,'prms'))
    {
      if(inherits(fixed_params,"mst_enorm"))
      {
        m = x$inputs$method
        x$inputs = x$mst_inputs
        x$inputs$method = m
        x$est$b = x$mst_est$b
        x$est$a = x$mst_est$a
      }
      
      if (fixed_params$inputs$method!="CML")
        message("Posterior means are taken as values for fixed parameters")
      
      fixed_params = fixed_params$inputs$ssIS
      fixed_params$b = if.else(fixed_params$inputs$method=="CML", fixed_params$est$b, colMeans(fixed_params$est$b))
      
    } else
    {
      # transform the fixed params to the b parametrization dexter uses internally
      fixed_params = transform.df.parms(fixed_params, out.format = 'b') 
    }  
    
    # avoid factor warnings and reduce
    fixed_params$item_id = factor(as.character(fixed_params$item_id), levels=levels(ss$design$item_id))
    fixed_params = filter(fixed_params,!is.na(.data$item_id)) 
    
    # check for missing categories in fixed_params
    missing_cat = ss$ssIS |> 
      semi_join(fixed_params, by='item_id') |>
      left_join(fixed_params, by=c('item_id','item_score')) |>
      filter(is.na(.data$b)) 
    
    if(nrow(missing_cat) > 0)
    {
      cat(paste('Some score categories are fixed while some are not, for the same item.',
                'Dexter does not know how to deal with that.\nThe following score categories are missing:\n'))
      missing_cat |> 
        select('item_id', 'item_score') |>
        arrange(.data$item_id, .data$item_score) |>
        as.data.frame() |>
        print()
      stop('missing score categories for fixed items, see output')
    }
  }
  
  if (method=="CML"){
    result = calibrate_CML(ss=ss, fixed_params=fixed_params)
  } else if(method=="Bayes")
  {
    result = calibrate_Bayes(ss, nIter=nDraws, fixed_params=fixed_params)
  } 
  
  est_mle = theta_wmle_c(b=matrix(if.else(method=="Bayes",colMeans(result$b), result$b),ncol=1),
                     a=ss$ssIS$item_score,
                     first=ss$design$first-1L, last = ss$design$last-1L,bk_nit = ss$booklet$nit, WLE=FALSE, n_cores=1L)

  mle = tibble(booklet_id=ss$booklet$booklet_id[est_mle$booklet], 
         booklet_score=drop(est_mle$booklet_score),
         theta=drop(est_mle$theta)) |>
    filter(is.finite(.data$theta))
  
  ss$method = method
  output = list(est=result, inputs=ss,abl_tables = list(mle = mle))
  
  class(output) = append('enorm', class(output)) 
  output
}

print.prms = function(x,...) print.enorm(x,...)

print.enorm = function(x, ...){
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


coef.prms = function(object, hpd = 0.95, what=c('items','var','posterior'), ...)
{
  coef.enorm(object, hpd = hpd, what=what, ...)
}

#' extract enorm item parameters
#' 
#' @param object an enorm parameters object, generated by the function \code{\link{fit_enorm}}
#' @param hpd width of Bayesian highest posterior density interval around mean_beta, 
#'  value must be between 0 and 1, default is 0.95 
#' @param what which coefficients to return. Defaults to \code{items} (the item parameters). Can also be \code{var} for the 
#' variance-covariance matrix (CML only) or \code{posterior} for all draws of the item parameters (Bayes only)  
#' @param ... further arguments to coef are ignored
#'  
#' @return 
#' Depends on the calibration method and the value of 'what'. For \code{what="items"}: 
#' 
#' \describe{
#' \item{CML calibration}{a data.frame with columns: item_id, item_score, beta, SE_beta}
#' \item{Bayesian calibration}{a data.frame with columns: item_id, item_score, mean_beta, SD_beta, <hpd_b_left>, <hpd_b_right>}
#' }
#' 
#' If \code{what="var"} or \code{what="posterior"} then  a matrix is returned with the variance-covariance matrix or the posterior draws 
#' respectively.
#' 
#' @details
#' 
#' The parametrisation of IRT models is far from uniform and depends on the author. Dexter uses the following parametrisation for the 
#' extended Nominal Response Model (NRM):
#' 
#' \deqn{
#' P(X=a_j|\beta,\theta) = \frac{\exp\left(a_j\theta-\sum_{g=1}^{j}\beta_g(a_g-a_{g-1})\right)}{1+\sum_h \exp\left(a_h\theta-\sum_{g=1}^{h}\beta_g(a_g-a_{g-1})\right)}
#' } 
#' 
#' where \eqn{a_j} is a shorthand for the integer score belonging to the j-th category of an item. 
#' 
#' For dichotomous items with \eqn{a_1=1} (i.e. the only possible scores are 0 and 1)
#' this formula simplifies to the standard Rasch model: \eqn{P(x=1|\beta,\theta)=\frac{\exp(\theta-\beta)}{1+\exp(\theta-\beta)}}. For polytomous items, 
#' when all scores are equal to the categories (i.e. \eqn{a_j=j} for all \eqn{j}) 
#' the NRM is equal to the Partial Credit Model, although with a different parametrisation than is commonly used. 
#' For dichotomous items and for all polytomous items where \eqn{a_j-a_{j-1}} is constant, the formulation is equal to the OPLM.
#' 
#' 
coef.enorm = function(object, hpd = 0.95, what=c('items','var','posterior'), ...)
{
  x = object
  what = match.arg(what)
  
  if(inherits(x,"mst_enorm"))
  {
    m = x$inputs$method
    x$inputs = x$mst_inputs
    x$inputs$method = m
  }
  
  if(what=='items')
  {
    if (x$inputs$method=="CML")
    {
      atab=data.frame(item_id=x$inputs$ssIS$item_id,
                      item_score=x$inputs$ssIS$item_score,
                      beta=x$est$beta,
                      SE_beta=sqrt(diag(x$est$acov.beta)),stringsAsFactors=FALSE)
    } else
    {
      if(hpd <= 0 ||  hpd >= 1)
        stop('hpd must be between 0 and 1')
      
      hh = t(apply(x$est$beta,2,hpdens, conf=hpd))
      atab=data.frame(item_id = x$inputs$ssIS$item_id,
                      a = x$inputs$ssIS$item_score,
                      mb = colMeans(x$est$beta),
                      sdb = apply(x$est$beta, 2, sd),
                      hpdl = hh[,1], hpdr=hh[,2],stringsAsFactors=FALSE)
      colnames(atab)=c("item_id" ,"item_score", "mean_beta", "SD_beta", 
                       sprintf("%i_hpd_b_left", round(100 * hpd)),
                       sprintf("%i_hpd_b_right", round(100 * hpd)))
    }
    rownames(atab) = NULL
    return(df_format(atab))
  } else if(what=='var')
  {
    if(x$inputs$method!="CML")
      stop('Variance-covariance matrix is only available for CML estimation')
    m = x$est$acov.beta
    colnames(m) = rownames(m) = paste(x$inputs$ssIS$item_id, x$inputs$ssIS$item_score)
    return(m)
  } else if(what=='posterior')
  {
    if(x$inputs$method!="Bayes")
      stop('The posterior of item parameters is only available for Bayesian estimation')
    m=x$est$beta
    colnames(m) = paste(x$inputs$ssIS$item_id, x$inputs$ssIS$item_score)
    return(m)
  }
}


# returns log likelihood, or vector of log likelihoods if bayes
logL = function(parms, mean_gibbs=FALSE)
{
  a = parms$inputs$ssIS$item_score
  b = parms$est$b
  if(is.matrix(b) && mean_gibbs)
    b=colMeans(b)
  
  scoretab = split(parms$inputs$scoretab,parms$inputs$scoretab$booklet_id,drop=FALSE)
  design = split(parms$inputs$design, parms$inputs$design$booklet_id, drop=TRUE)
  
  llb = function(b)
  {
    ll_rm = sum(parms$inputs$ssIS$sufI * log(b))
    
    lgRM = sapply(scoretab,
                  function(stb)
                  {
                    d = design[[stb$booklet_id[1]]]
                    # impossible scores get a result of 0
                    # since they always have N=0 by definition
                    # it's easier to filter them out
                    w = stb$N>0
                    # c function, 0 index
                    log(elsymC(b, a, d$first-1L, d$last-1L)[w]) * stb$N[w]
                  }) |>
      unlist() |>
      sum()
    
    ll_rm-lgRM
  }
  if(is.matrix(b))
  {
    apply(b,1, llb)
  } else
  {
    llb(b)
  }
}

logLik.enorm = function(object,...)
{
  ll = logL(object)
  
  attr(ll, "df") = nrow(object$inputs$ssIS)-1L
  class(ll) = "logLik"
  ll
}




# Calibration -------------------------------------------------------------

lbnm = function(max_score)
{
  out = matrix(0,max_score+1,max_score+1)
  for(s in 2:max_score)
  {
    out[2:(s+1),s+1] = lchoose(s,1:s)
  }
  out + t(out)
}


implicit_update_b = function(b,indx,EsufI,sufI, n_rsp, item_fixed)
{
  max_dif = 0
  for(i in seq_along(indx)) if(!item_fixed[i])
  {
    x = indx[[i]]

    e0 = c(max(n_rsp[i] - sum(EsufI[x]),1), EsufI[x])
    i0 = c(n_rsp[i]-sum(sufI[x]),sufI[x])
    
    b_item = b[x]
    new_b_item = c(1,b_item) * i0/e0
    b[x] = new_b_item[-1]/new_b_item[1]
    
    max_dif = max(max_dif, abs(e0-i0))
  }
  list(max_dif=max_dif,b=b)
}


# include the 0 category in determining convergence
EsufI_dif = function(indx,EsufI,sufI, n_rsp, item_fixed)
{
  max_dif = 0
  for(i in seq_along(indx)) if(!item_fixed[i])
  {
    x = indx[[i]]
    
    e0 = c(max(n_rsp[i] - sum(EsufI[x]),1), EsufI[x])
    i0 = c(n_rsp[i]-sum(sufI[x]),sufI[x])

    max_dif = max(max_dif, abs(e0-i0))
  }
  max_dif
}


calibrate_CML = function(ss,fixed_params=NULL) 
{
  max_ie_iter = 1000
  max_nr_iter = 30
  use_mean = FALSE
  n_cores = get_ncores(min(nrow(ss$booklet), 64))
  
  pb = get_prog_bar()
  on.exit({pb$close()})

  # preparation
  lbinom = NULL
  
  # THESE ARE C INDEXED!
  # one vector for all booklets
  bkfirst = as.integer(ss$design$first - 1L)
  bklast = as.integer(ss$design$last - 1L)
  
  nbk = nrow(ss$booklet)
  nit = nrow(ss$ssI)
  
  sufI = ss$ssIS$sufI
  a = ss$ssIS$item_score
  
  indx = mapply(ss$ssI$first, ss$ssI$last, FUN=':',SIMPLIFY=FALSE)
  
  
  if(is.null(fixed_params))
  {
    b = rep(1,length(a))
    nn = sum(ss$ssI$n_rsp)
    item_fixed = logical(nrow(ss$ssI))
    fixed_b = NULL
  } else
  {
    b = fixed_params |>
      right_join(ss$ssIS, by=c('item_id','item_score')) |>
      arrange(.data$item_id,.data$item_score) |>
      pull(.data$b)
    
    if(!anyNA(b)) stop('nothing to calibrate, all parameters are fixed')
    
    fixed_b = !is.na(b)
    b = coalesce(b,1)
    
    item_fixed = (ss$ssI$item_id %in% fixed_params$item_id)
    nn = sum(ss$ssI$n_rsp[item_fixed==0])
  }
  
  
  if(use_mean)
    lbinom = lbnm(max(ss$booklet$max_score)+3)
  
  ## Implicit Equations  ###
  converged=FALSE
  
  ie_iter = 0
  while (!converged && ie_iter <= max_ie_iter)
  {
    ie_iter = ie_iter+1
    
    if(use_mean)
      EsufI = Expect_binom(lbinom, b, a, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
    else
      EsufI = Expect(b, a, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
    
    
    res = implicit_update_b(b,indx,EsufI,sufI, ss$ssI$n_rsp,item_fixed)

    converged = (res$max_dif/nn) < 1e-04
    
    if(is.na(converged))
    {
      if(use_mean)
        stop('problem cannot be computed')
      
      use_mean = TRUE
      converged = FALSE
      lbinom = lbnm(max(ss$booklet$max_score)+3)
      
      next
    } 
    b = res$b
    pb$tick()
  }
  if (!converged) warning(paste('Implicit Equations not Converged in',ie_iter,"iterations"))
    
  ### identification ###
  if(is.null(fixed_b))
  {
    ref_cat = 1
    b = b/(b[ref_cat]^(a/a[ref_cat]))
  }
    
  ###  NR  ###
  H = matrix(0,length(a),length(a))
  converged = FALSE
  nr_iter=0

  while (!converged && nr_iter<max_nr_iter)
  {
    nr_iter=nr_iter+1
    
    if(use_mean)
      Hess_binom(lbinom,b, a, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, n_cores, EsufI, H)
    else
      Hess(b, a, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, n_cores, EsufI, H)
    
    if(is.null(fixed_b))
    {
      # identify
      H[ref_cat,]=0
      H[,ref_cat]=0
      H[ref_cat,ref_cat]=1
      EsufI[ref_cat] = sufI[ref_cat] 
    } else
    {
      H[fixed_b,]=0
      H[,fixed_b]=0
      diag(H)[fixed_b]=1
      EsufI[fixed_b]=sufI[fixed_b]
    }
    
    # take 0 categrory into account in determining convergence
    converged = EsufI_dif(indx, EsufI, sufI, ss$ssI$n_rsp, item_fixed)/nn <1e-10
      
    nb = try(b*exp(solve(H,sufI-EsufI)))
    
    if(inherits(nb,'try-error') || is.na(converged))
    {
      if(use_mean)
        stop('problem cannot be computed')
      
      use_mean = TRUE
      converged = FALSE
      lbinom = lbnm(max(ss$booklet$max_score)+3)
      next
    }
    
    if(is.null(fixed_b))
    {
      b[-ref_cat] = nb[-ref_cat]
    } else
    {
      b[!fixed_b] = nb[!fixed_b]
    }
    pb$tick()
  }
  if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  
  report = toOPLM(ss$ssIS$item_score, b, ss$ssI$first, ss$ssI$last, H=H, fixed_b=fixed_b)
  b = report$b_renorm
  
  lx = mapply(function(bdes, bsct)
      {
        g = elsymC(b, a, bdes$first-1L, bdes$last-1L) # c function, zero index
        tibble(booklet_id = bdes$booklet_id[1],
               booklet_score = bsct$booklet_score,
               lambda = ifelse(bsct$N>0, bsct$N/g, NA_real_)) 
      },
      split(ss$design, ss$design$booklet_id),
      split(ss$scoretab, ss$scoretab$booklet_id),
      SIMPLIFY=FALSE, USE.NAMES=FALSE) |>
      bind_rows()
  
  return(list(b=b, H=H, beta=report$beta, acov.beta=report$cov.beta, lambda=lx,
               n_iter=ie_iter+nr_iter, nr_iter = nr_iter, ie_iter=ie_iter))
  
}



calibrate_Bayes = function(ss,  nIter, fixed_params=NULL, 
                           from = 50L, step = 10L,
                           start_b=NULL)
{
  nchains = ncores = get_ncores(desired = min(ceiling(nIter/4),128L), maintain_free = 1L)

  pb = get_prog_bar(nsteps = nIter)
  on.exit({pb$close()})
  
  # decide on starting values
  if(!is.null(start_b))
  {
    b = matrix(rep(start_b,nchains), ncol=nchains)
  } else
  {
    cml = calibrate_CML(ss, fixed_params)
    
    if(nchains==1)
    {
      b = matrix(cml$b,ncol=1)
    } else
    {
      # different but slightly underdispersed start locations
      sample_beta = rmvnorm(nchains,cml$beta, cml$acov.beta/2)
      
      b = apply(sample_beta,1,function(beta){ beta2b(ss$ssIS$item_score,beta,ss$ssI$first, ss$ssI$last)})
      
    }
  }
  
  if(is.null(fixed_params))
  {
    item_fixed = integer(nrow(ss$ssI))
    fixed_b=NULL
  } else
  {
    item_fixed = as.integer(ss$ssI$item_id %in% fixed_params$item_id)
    fixed_b = fixed_params |>
      right_join(ss$ssIS, by=c('item_id','item_score')) |>
      arrange(.data$item_id,.data$item_score) |>
      pull(.data$b)
    
    # prevent rounding/etc errors in strating values
    b[!is.na(fixed_b),] = fixed_b[!is.na(fixed_b)]
  }

  prior_eta = 0.5
  prior_rho = 0.5
  prior_nu = 0.1
  
  # bookkeeping: make some extra counts and indexes for C function
  design = ss$design |>
    mutate(bn = dense_rank(.data$booklet_id) - 1L) |>
    group_by(.data$bn) |>
    mutate(inr = dense_rank(.data$first) - 1L) |>
    ungroup() |>
    arrange(.data$first, .data$bn)
  
  bi = design$inr
  ib = design$bn
  
  nbi = design |> 
    count(.data$first) |> 
    arrange(.data$first) |> 
    pull(.data$n)
  
  design = arrange(design,.data$bn)
  
  bfirst = as.integer(design$first -1L)
  blast = as.integer(design$last -1L)
  
  #end bookkeeping
  
  out = calibrate_Bayes_chains(ss$ssIS$item_score, ss$ssI$first-1L, ss$ssI$last-1L,
                                 ib, bi, nbi, ss$booklet$nit, bfirst, blast, ss$booklet$max_score, ss$booklet$M,
                                 ss$ssIS$sufI, ss$ssI$n0, ss$scoretab$N, b, item_fixed, 
                                 from, step, 
                                 as.integer(nIter), pb$cpp_prog_init(), ncores,
                                 prior_eta, prior_rho, prior_nu)


  
  
  report = toOPLM(ss$ssIS$item_score, out$b, ss$ssI$first, ss$ssI$last, H=NULL,method='Bayes',fixed_b=fixed_b)

  colnames(out$lambda) = paste(ss$scoretab$booklet_id, ss$scoretab$booklet_score, sep='-')
  
  return(list(a=ss$ssIS$item_score, b=report$b_renorm, 
              lambda=out$lambda, beta=report$beta,gibbs_b=out$b,
              chain_start=out$chain_start,
              priors=list(eta=prior_eta, rho=prior_rho, nu=prior_nu))) 
}






