
#' Estimate the Interaction and the Rasch model
#'
#' Estimate the parameters of the Interaction model and the Rasch model
#'
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @return An object of class \code{inter} holding results
#' for the Rasch model and the interaction model.
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function therefore fits the
#' Rasch model and the interaction model on complete data. 
#' This typically consist of responses to items in one booklet but can also consist of
#' the intersection (common items) in two or more booklets. If the intersection is empty
#' (no common items for all persons), the function will exit with an error message.
#'
#' @seealso \code{\link{plot.inter}}, \code{\link{fit_domains}}
#'
#' @examples
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules, ":memory:")
#' add_booklet(db, verbAggrData, "agg")
#'
#' m = fit_inter(db, booklet_id=='agg')
#' plot(m, "S1DoScold", show.observed=TRUE)
#'
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#'
fit_inter = function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  check_dataSrc(dataSrc)
  fit_inter_(dataSrc, qtpredicate, env, regs=TRUE)
}

fit_inter_ = function(dataSrc, qtpredicate = NULL, env=NULL, regs=TRUE)
{
  
  respData = get_resp_data(dataSrc, qtpredicate, env = env, retain_person_id=FALSE) |>
	  intersection_rd()

  if(nrow(respData$x)==0) 
    stop('no responses to analyse')

  ss = get_sufStats_im(respData)

  

  est = calibrate_rim(ss, regs=regs)

  if(regs)
  {
    est$itrRM = rowsum(est$ctrRM * ss$ssIS$item_score, ss$ssIS$item_id, reorder=FALSE)
    est$itrIM = rowsum(est$ctrIM * ss$ssIS$item_score, ss$ssIS$item_id, reorder=FALSE)
  }
  
  ss$design=respData$design
  output = list(est = est, inputs = ss)
  class(output) = append("inter", class(output))
  output
}


#' Estimate the Rasch and the Interaction model per domain
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param item_property The item property defining the
#' domains (subtests)
#' @return An object of class \code{imp} holding results
#' for the Rasch model and the interaction model.
#' @details 
#' We have generalised the interaction model for items having more than two (potentially, a largish number) 
#' of response categories. This function represents scores on subtests as 
#' super-items and analyses these as normal items.
#'
#' @seealso \code{\link{plot.inter}}, \code{\link{fit_inter}}, \code{\link{add_item_properties}}
#'
#' @examples
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules, ":memory:")
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' mSit = fit_domains(db, item_property= "situation")
#' plot(mSit)
#'
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#'
fit_domains = function(dataSrc, item_property, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()

  get_resp_data(dataSrc, qtpredicate, extra_columns = item_property, env = env, retain_person_id=FALSE) |>
	  intersection_rd() |>
    polytomize_rd(item_property, protect_x = !is_db(dataSrc)) |>
    fit_inter_()

}



print.inter = function(x, ...){
  res = paste0('Parameters for the Rasch and Interaction Model', 
               '\n\n# use plot() for plotting the Rasch and Interaction Model or coef() for retreiving the parameters\n')
  cat(res)
  invisible(res)
}

#' Extract interaction model parameters
#' 
#' @param object an object returend by the function \code{\link{fit_inter}}
#' @param what which coefficients to return. Defaults to \code{items} (the item parameters), can also be \code{scoreprob}
#' for the probability of each item score per booklet score.
#' @param ... further arguments to coef are ignored
#' 
coef.inter = function(object, what=c("items","scoreprob"), ...) 
{
  x = object
  what = match.arg(what)
  if(what == 'items')
  {
    first = x$inputs$ssI$first
    last  = x$inputs$ssI$last
    report_RM = toOPLM(x$inputs$ssIS$item_score, x$est$bRM, first, last)
    report_IM = toOPLM(x$inputs$ssIS$item_score, x$est$bIM, first, last)
    
    IS = tibble(item_id = x$inputs$ssIS$item_id, item_score = x$inputs$ssIS$item_score,
                beta_rasch = as.vector(report_RM$beta), beta_IM = as.vector(report_IM$beta))
    I = tibble(item_id = x$inputs$ssI$item_id, sigma = log(x$est$cIM), SE_sigma= x$est$se.sigma, fit_IM=x$est$fit.stats)
    
    inner_join(IS,I,by='item_id') |> 
  	  arrange(.data$item_id, .data$item_score) |> 
      mutate(item_id=as.character(.data$item_id)) |>
  	  df_format()
  } else
  {
    z = x$est$ctrIM

    tibble(item_id = rep(as.character(x$inputs$ssIS$item_id),ncol(z)), 
           item_score = rep(x$inputs$ssIS$item_score,ncol(z)),
           booklet_score = rep(x$inputs$scoretab$booklet_score,each=nrow(z)),
           p = as.double(z)) |>
      df_format()
  }
}



calibrate_rim = function(ss, regs=FALSE) {
  
  ssI = ss$ssI
  ssIS = ss$ssIS
  nit=nrow(ssI)
  scoretab = ss$scoretab$N

  a = ss$ssIS$item_score
  m = sum(scoretab)

  b = rep(1,nrow(ssIS))
  ic = rep(1,nit)
  var.ic = double(nit)
  HIM = vector("list", nit)
  
  
  first0 = ssI$first-1L
  last0 = ssI$last-1L
  
  ps = possible_scores(a, ssI$first, ssI$last)
  # see actual scores, do: (1:length(ps)-1L)[as.logical(ps)]
  
  converged=2
  scale=2
  itr_rasch=0L
  while(converged>0.01)
  {
    converged=-1
    pi_mat = ittotmatC(b,a,ic,first0,last0, ps)
    
    for (i in 1:nit)
    {
      w = ssI$first[i]:ssI$last[i]
      pi = pi_mat[w,,drop=FALSE]
      E = ssIS$sufI[w]- pi%*%scoretab
      H = -pi %*% tcrossprod(diag(scoretab), pi) 
      diag(H) = pi%*%scoretab + diag(H)
        
      # NR update for parameters of item i
      update = solve(H*scale,E)
      b[w]=b[w]*exp(update)
      converged=pmax(converged,max(abs(E))/m)
    }
    if (converged<1) scale=1
    itr_rasch=itr_rasch+1L
  }
  
  bRM = b

  ## IM
  converged=2
  scale=2
  first_iter = TRUE
  itr_IM=0L
  while(converged>0.001)
  {
    converged=-1
    pi_mat = ittotmatC(b,a,ic,first0,last0, ps)
    
    if(first_iter)
    {
      # save final pi_mat for Rasch
      first_iter = FALSE
      ctrRM = pi_mat
    }
    for (i in 1:nit)
    {
      # gradient and hessian for thresholds of item i
      w = ssI$first[i]:ssI$last[i]
      pi = pi_mat[w,,drop=FALSE]
      E = ssIS$sufI[w] - pi %*% scoretab
      H = -pi %*% tcrossprod(diag(scoretab), pi)
      diag(H) = pi %*% scoretab + diag(H)
        
      # gradient and hessian for interaction parameter
      ncol_pi = ncol(pi)
      nrow_pi = nrow(pi)
      s_range = 0:(ncol_pi-1)
        
      E = c(E, ssI$sufC[i])
      H = cbind(H, rep.int(0, nrow(H)))
      H = rbind(H, rep.int(0, ncol(H)))
      k = 1
      e0 = 0; e1 = 0
      f = matrix(0, nrow_pi, ncol_pi)
      h = 0
      for (j in w)
      {
        E[length(E)] = E[length(E)]-a[j]*sum(s_range*scoretab*pi[k,])
        e0=e0+a[j]*pi[k,]
        e1=e1+a[j]^2*pi[k,]
        f[k,]=a[j]*s_range*pi[k,]
        h=h+f[k,]
        k=k+1
      }
      H[nrow(H),nrow(H)]=sum(s_range^2*(e1-e0^2)*scoretab)
      for (k in 1:nrow(f))
      {
        H[k,nrow(H)]=sum((f[k,]-pi[k,]*h)*scoretab)
        H[nrow(H),k]=H[k,nrow(H)]
      }
        
        
      # NR update for parameters of item i
      update = solve(H*scale,E)
      b[w] = b[w]*exp(update[-length(update)])
      ic[i] = ic[i]*exp(update[length(update)])
      HIM[[i]] = H
      var.ic[i] = solve(H)[nrow(H),nrow(H)]
      converged = pmax(converged,max(abs(E))/m)
      
    }
    if (converged<1) scale=1
    itr_IM=itr_IM+1L
  }
  
  sigma = log(ic)
  sigma = sigma - mean(sigma)
  ic = exp(sigma)
  se.sigma = sqrt(var.ic)
  fit.stats = sigma/se.sigma
  
  cIM_score = double(length(b))
  for(i in 1:nrow(ssI)) 
    for(j in ssI$first[i]:ssI$last[i]) 
      cIM_score[j] = ic[i]
  
  
  out = list(bRM=bRM,cRM=rep(1,nit),bIM=b,cIM=ic,cIM_score=cIM_score, se.sigma=se.sigma,HIM=HIM, fit.stats=fit.stats, 
             possible_scores = (1:length(ps)-1L)[as.logical(ps)],
             itr_rasch=itr_rasch,itr_IM=itr_IM)
  if(regs)
  {
    out$ctrRM = ctrRM
    out$ctrIM = ittotmatC(b,a,ic,first0,last0, ps) 
  }

  out
}






#' Simulation from the interaction model
#'
#' Simulate item scores conditional on test scores using the interaction model
#'
#' @param m an object produced by function \code{fit_inter}
#' @param scores vector of test scores
#' 
#' @return
#' a matrix with item scores, one column per item and one row per test score. Row order
#' equal to scores
#' 
r_score_IM = function(m, scores)
{

  if(inherits(m,'data.frame'))
  {
    stop('input `m` must be of class "inter"')
    # this does not yet work
    # if('beta_IM' %in% colnames(m) && !'beta' %in% colnames(m))
    #   m$beta = m$beta_IM
    # m = arrange(m,.data$item_id, .data$item_score)
    # prms = simplify_parms(m)
    # 
    # a = prms$a
    # bIM = prms$b
    # first = prms$items$first
    # last = prms$items$last
    # cIM = m$sigma
    
  } else if(inherits(m,'inter'))
  {
    first = m$inputs$ssI$first
    last = m$inputs$ssI$last
    a = m$inputs$ssIS$item_score
    bIM = m$est$bIM
    cIM = m$est$cIM
  } 
  else stop('input `m` must be of class "inter"')
 
  maxs = sum(a[last])
  
  if(any(scores>maxs))
    stop('scores may not be larger than the maximum score on the test')
  
  if(any(scores<0))
    stop('all scores must be positive')
  
  
  scoretab = score_tab_single(scores, maxs)
  s = sampleIMC(bIM,cIM,a,as.integer(first-1L), as.integer(last-1L), scoretab)
  
  if(scoretab[1]>0)
    s[1:scoretab[1],] = 0
  
  if(scoretab[maxs+1]>0)
    s[(nrow(s)-scoretab[maxs+1]+1):nrow(s),] = matrix(a[last], nrow=scoretab[maxs+1], ncol = ncol(s), byrow=TRUE)
  
  colnames(s) = m$inputs$ssI$item_id
  
  if(is.unsorted(scores))
    s = s[order(order(scores)),,drop=FALSE]

  s
}