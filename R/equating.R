
# to do: example    
# to do: parm design, describe what it is for and only then what it should look like
# parm design migfht also accept a vector of item id's?    
#' The probability to pass on a reference test given a score on a new booklet
#' 
#' Given response data that form a connected design,
#' compute the probability to pass on the reference set conditional on each score on one or more target tests.
#' 
#' Note that this function is computationally intensive and can take a long time to run, especially when computing the
#' probability to pass for multiple target booklets. 
#'
#' @param dataSrc Data source: a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score and booklet_id
#' @param ref_items vector with id's of items in the reference set, they must all occur in dataSrc
#' @param pass_fail pass-fail score on the reference set, the lowest score with which one passes
#' @param predicate An optional expression to subset data in dataSrc, if NULL all data is used
#' @param design A data.frame with columns booklet_id (if multiple booklets) and item_id defining the target test booklet(s), 
#' if NULL (default) this will be derived from the dataSrc and the probability to pass will be computed 
#' for each test score for each booklet in your data.
#' @param similar_groups When TRUE it is assumed that candidates taking the reference test and the target tests are similar 
#' in ability. If they are similar, the score distribution on the target test is estimated more precisely. 
#' Works only when design contains booklet_id's.
#' 
#' @return An object of type p2pass. Use \code{coef()} to extract the 
#' probablity to pass for each booklet and score. Use \code{plot()} to plot 
#' the probabilities, sensitivity and specificity. 
#'         
#' @details
#' For any possible score k, we use a Gibbs sampler to calculate 
#' \ifelse{html}{
#'   \out{<p><i>P(Y<sub>+</sub> &#8805; c | X<sub>+</sub>=k, <span style="font-weight:bold;">x</span> ) = 
#'              <span style="font-size:150\%;">&#8747;</span><sub><span style="font-weight:bold;">b</span>,&theta;</sub> 
#'              P(Y<sub>+</sub> &#8805; c | &theta;, <span style="font-weight:bold;">b</span>)&#8729;&#402;(<span style="font-weight:bold;">b</span> | <span style="font-weight:bold;">x</span>)&#8729;d<span style="font-weight:bold;">b</span>&theta;
#'              </i></p>}
#' }{\deqn{P(Y_+ \geq c|X_+ = k,\mathbf{x})
#'      = \int_{\mathbf{b},\theta} 
#'             P(Y_+ \geq c|\theta,\mathbf{b}) 
#'             f(\theta|X_+ = k) f(\mathbf{b}|\mathbf{x})
#'      d\mathbf{b},\theta}} 
#' 
#' where:
#' \describe{
#'    \item{\eqn{\theta}}{is student ability}
#'    \item{\eqn{b}}{are the item parameters}
#'    \item{\ifelse{html}{\out{<i>Y<sub>+</sub></i>}}{\eqn{Y_+}}}{is the score on the reference test}
#'    \item{\eqn{c}}{is an established pass_fail score on the reference test}
#'    \item{\ifelse{html}{\out{<i>X<sub>+</sub></i>}}{\eqn{X_+}}}{is the score on the booklet}
#'    \item{\eqn{x}}{are the observed data}
#'  }
#' This probability can be used to establish a pass-fail score for the new booklet.     
#'    
probability_to_pass = function(dataSrc, ref_items, pass_fail, design = NULL, predicate = NULL, similar_groups = TRUE)
{
  check_dataSrc(dataSrc)
  check_num(pass_fail, 'integer', .min=0)
  check_df(design, 'item_id', nullable=TRUE)
  
  
  # optionally these can become extra parameters
  est_method = "Bayes"
  nIterations = 1000 ## number over which we average
  include_theta = TRUE
  ##
 
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env)
  
  ref_items = tibble(item_id = ffactor(unique(ref_items), levels=levels(respData$design$item_id)))
  
  if(anyNA(ref_items$item_id))
    stop('One or more of your reference items does not occur in your dataSrc')

  if(is.null(design))
  {
    design = respData$design
  } else
  {
    design$item_id = ffactor(as.character(design$item_id), levels = levels(respData$design$item_id))
    if(anyNA(design$item_id))
      stop('One or more items in design does not occur in your dataSrc ')
    
    if(!('booklet_id' %in% colnames(design))) 
      design$booklet_id = 'target items'
  }
  
  # Use all available data to estimate the item parameters
  f = fit_enorm_(respData, method=est_method, nIterations = 5*nIterations)
  
  # Now reduce data 
  items_rel = union(ref_items$item_id, design$item_id)
  bkl_rel = respData$design %>% 
    mutate(item_id = as.character(.data$item_id)) %>%
    inner_join(tibble(item_id=items_rel),by='item_id') %>% 
    distinct(.data$booklet_id)
  
  respData = semi_join(respData, bkl_rel, by='booklet_id') #on booklet, so no recompute sumscores
  
  ref_ssI = f$inputs$ssI %>% 
    semi_join(ref_items, by = 'item_id') %>% 
    arrange(.data$first)
  
  ref_first = pull(ref_ssI, 'first')
  ref_last = pull(ref_ssI, 'last')
  
  # Get mean and sd of ability in sample
  #ab = ability(respData, f, method="MLE") 
  #new_mu = mean(ab$theta[is.finite(ab$theta)])
  #new_sigma = sd(ab$theta[is.finite(ab$theta)])
  pv = plausible_values(respData, f)
  new_mu = mean(pv$PV1)
  new_sigma = sd(pv$PV1)
  
  
  bkl = unique(design$booklet_id)
  bk_results = lapply(bkl, function(booklet)
  {
    dsg = filter(design, .data$booklet_id == booklet)
    ssI = semi_join(f$inputs$ssI, dsg, by = 'item_id') %>% arrange(.data$first)
    a = pull(f$inputs$ssIS, 'item_score')
    b = f$est$b
    
    new_first = pull(ssI,'first')
    new_last = pull(ssI,'last')
    scores = 0:sum(a[new_last])
    # to do: spurious equal booklet names??
    if (similar_groups || booklet=='target booklet')
    {
      p_new = plausible_scores(respData, parms=f, items = design$item_id) %>%
        count(.data$PS1) %>%
        full_join(tibble(PS1 = as.integer(0:sum(a[ssI$last]))), by='PS1') %>%
        mutate(n = coalesce(n,0L)/sum(n, na.rm=TRUE)) %>%
        pull(.data$n)
    }else
    {
        rel_bkl = respData$design %>% filter(.data$booklet_id==booklet)
        rel_resp = semi_join(respData, rel_bkl, by='booklet_id')
        p_new = plausible_scores(rel_resp, parms=f, items = pull(dsg, 'item_id')) %>%
          group_by(.data$PS1) %>%
          summarise(n = as.integer(n())) %>%
          ungroup() %>%
          full_join(tibble(PS1 = as.integer(0:sum(a[ssI$last]))), by='PS1') %>%
          mutate(n = coalesce(n,0L)/sum(n, na.rm=TRUE)) %>%
          pull(.data$n)
    }
    
    iter_set = seq(5, 5*nIterations, by=5)
    if (est_method == "Bayes")
    {
      probs = matrix(0, sum(a[ssI$last])+1, length(iter_set))
      ref_range = (pass_fail:sum(a[ref_ssI$last])) + 1
      eq_score = rep(NA, length(iter_set))
      
      tel = 1
      pb = txtProgressBar(min=0, max=length(iter_set))
      for (iter in iter_set)
      {
        if (include_theta)
        {
          ref_theta = theta_MLE(b[iter,],a, ref_first, ref_last)$theta[pass_fail+1]
          suppressWarnings({eq_score[tel] = min(which(theta_MLE(b[iter,],a, new_first, new_last)$theta >= ref_theta))-1 })
          if(is.infinite(eq_score[tel])) browser()
        }
        spv = pv_recycle(b[iter,], a, new_first, new_last, scores, npv=1, mu=new_mu, sigma=new_sigma)
        #spv = recycle_pv(b[iter,], a, new_first, new_last, npv=1, mu=new_mu, sigma=new_sigma)
        #prf = apply(spv,1, function(x)pscore(x,b[iter,],a,ref_first,ref_last)) 
        prf = pscore(spv, b[iter,],a,ref_first,ref_last)
        probs[,tel] = apply(prf, 2, function(x) sum(x[ref_range]))
        tel=tel+1
        setTxtProgressBar(pb, value=tel)
      }
      close(pb)
    } else
    {
      probs = matrix(0,sum(a[new_last])+1,nIterations)
      ref_range = (pass_fail:sum(a[ref_last]))+1
      ref_theta = theta_MLE(b,a,ref_first,ref_last)$theta[pass_fail+1]
      eq_score = min(which(theta_MLE(b,a,new_first,new_last)$theta >= ref_theta))-1 
      
      spv = pv_recycle(b, a, new_first, new_last, scores, npv=nIterations, mu=new_mu, sigma=new_sigma)
      for (iter in 1:nIterations)
      {
        prf = pscore(spv[,iter], b, a, ref_first, ref_last)
        probs[,iter] = apply(prf, 2, function(x) sum(x[ref_range]))
      }
    }

    
    p_pass_given_new = rowMeans(probs) 
    ## additional
    tp = rev(cumsum(rev(p_pass_given_new*p_new)))/rev(cumsum(rev(p_new)))
    tp_rate = rev(cumsum(rev(p_pass_given_new*p_new)))/tp[1]
    tn = rev(cumsum(rev((1-p_pass_given_new)*p_new)))/rev(cumsum(rev(p_new)))
    fp_rate =  rev(cumsum(rev((1-p_pass_given_new)*p_new)))/tn[1]
    
    # to do: this output format seems somewhat inefficient as 3 data.frames all share a column
    # might make it 1 df?
    
    list(prob_to_pass = 
           tibble(score_new = 0:sum(a[new_last]), probability_to_pass = p_pass_given_new, true_positive = tp),
         properties = 
           tibble(score_new = 0:sum(a[new_last]), sensitivity = tp_rate, specificity = 1-fp_rate),
         roc = 
           tibble(score_new = 0:sum(a[new_last]), true_positive_rate = tp_rate, false_positive_rate = fp_rate),
         pnew = p_new, 
         equated_score=eq_score)
  })
  names(bk_results) = bkl

  out = list(booklets = bk_results, est_method = est_method)
  class(out) = append('p2pass', class(out))
  out
}


print.p2pass = function(x,...)
{
  cat(paste('Equating information for', length(x$booklets),
              'booklets. Use `coef()` to extract the statistics or `plot()` to',
              'plot the probabilty to pass and sensitivity/specificity.\n'))
  invisible(x)
}


#' extract equating information
#' 
#' @param object an p2pass object, generated by \code{\link{probability_to_pass}}
#' @param ... further arguments are currently ignored
#' @return 
#' A data.frame with columns:
#' \describe{
#' \item{booklet_id}{id of the target booklet}
#' \item{score_new}{score on the target booklet}
#' \item{probability_to_pass}{probability to pass on the reference test given score_new}
#' \item{true_positive}{percentages that correctly passes}
#' \item{sensitivity}{The proportion of positives that are correctly identified as such}
#' \item{specificity}{The proportion of negatives that are correctly identified as such}
#' \item{proportion}{proportion in sample with score_new}}
coef.p2pass = function(object, ...)
{
  # this join might be better done in p2pass
  lapply(object$booklets, function(bk)
    {
      bk$prob_to_pass %>%
        inner_join(bk$properties, by='score_new') %>%
        add_column(probability = bk$pnew)
    }
  ) %>%
    bind_rows(.id = 'booklet_id') %>%
    as.data.frame()
}
# to do: do we also want the roc plot?
##########################################
#' A plot method for probability_to_pass
#'
#' Plot equating information from probability_to_pass 
#'
#'
#' @param x An object produced by function \code{\link{probability_to_pass}}
#' @param booklet_id vector of booklet_id's to plot, if NULL all booklets are plotted
#' @param what information to plot, 'equating', 'sens/spec' or 'both'
#' @param ... Any additional plotting parameters.
#' @method plot p2pass
#'
plot.p2pass = function(x, ..., booklet_id=NULL, what = c('both','equating','sens/spec'))
{
  user.args = list(...)
  what = match.arg(what)
  
  if(!is.null(booklet_id))
  {
    bk_results = x$booklets[booklet_id] # as list on purpose
  } else
  {
    bk_results = x$booklets
  }
  
  for(bkl in names(bk_results))
  {
    if(what != 'sens/spec')
    {
      #equating plot
      dflt = list(xlab = paste("score on", bkl), ylab = "Probability to pass", main = "Equating",
                 pch = 16, bty='l', col='grey40')
      ovr = list(x = bk_results[[bkl]]$prob_to_pass$score_new, 
                 y = bk_results[[bkl]]$prob_to_pass$probability_to_pass)
      do.call(plot, merge_arglists(user.args, dflt, ovr))
      if(x$est_method == 'Bayes')
      {
        eql = table(bk_results[[bkl]]$equated_score) / length(bk_results[[bkl]]$equated_score)
        xc = as.integer(names(eql))
        yc = as.vector(eql)
        segments(xc,0,xc,yc, col='red')
      }
      
      
    }
    if(what != 'equating')
    {
      # sens/spec plot
      dflt = list(xlab = paste("score on", bkl),main="Sensitivity/Specificity", ylab = "sens./spec.",
                  pch = 16, ylim=c(0,1), bty='l',col='#1a9641')
      ovr = list(x = bk_results[[bkl]]$properties$score_new, 
                 y = bk_results[[bkl]]$properties$sensitivity,
                 type="o")
      do.call(plot, merge_arglists(user.args, dflt, ovr))
      
      
      lnargs = dropNulls(merge_arglists(
        user.args[c('lty','lwd','col','pch')],
        list(col='#d7191c',pch = 16),
        list(type="o",
             x = bk_results[[bkl]]$properties$score_new,
             y = bk_results[[bkl]]$properties$specificity)))
        
      do.call(lines, lnargs)
    }
  }
  
  invisible(NULL)
}

