

# to~do: example    
  

#' The probability to pass on a reference test given a score on a new booklet
#' 
#' Given response data that form a connected design,
#' compute the probability to pass on the reference set conditional on each score on one or more target tests.
#' 
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms parameters returned from fit_enorm. If uncertainty about parameter estimation should be included
#' in the computations, use `method='Bayes'` and nDraws equal or larger than nDraws in probability_to_pass
#' @param ref_items vector with id's of items in the reference set, they must all occur in dataSrc
#' @param pass_fail pass-fail score on the reference set, the lowest score with which one passes
#' @param predicate An optional expression to subset data in dataSrc, if NULL all data is used
#' @param target_booklets The target test booklet(s). A data.frame with columns booklet_id (if multiple booklets) and item_id, 
#' if NULL (default) this will be derived from the dataSrc and the probability to pass will be computed 
#' for each test score for each booklet in your data.
#' @param nDraws The function uses an Markov-Chain Monte-Carlo method to calculate the probability to pass and this is the number of Monte-Carlo samples used. 
#' @return An object of type \code{p2pass}. Use \code{coef()} to extract the 
#' probablity to pass for each booklet and score. Use \code{plot()} to plot 
#' the probabilities, sensitivity and specificity or a ROC-curve. 
#'         
#' @details
#' Note that this function is computationally intensive and can take some time to run, especially when computing the
#' probability to pass for multiple target booklets. Further technical details can be found in a vignette.
#' 
#' @seealso The function used to plot the results: \code{\link{plot.p2pass}}

probability_to_pass = function(dataSrc, parms, ref_items, pass_fail, predicate = NULL, 
                                 target_booklets = NULL, nDraws = 1000)
{
  check_dataSrc(dataSrc)
  check_num(pass_fail, 'integer', .min=0, .length=1)
  check_df(target_booklets, 'item_id', nullable=TRUE)
  design=target_booklets

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
  n_booklets = n_distinct(design$booklet_id)
  pgb = prog_bar((nDraws+2) * n_booklets)
  on.exit({pgb$close()})

  ref_ssI = parms$inputs$ssI %>% 
    semi_join(ref_items, by = 'item_id') %>% 
    arrange(.data$first)
  
  ref_first = ref_ssI$first
  ref_last = ref_ssI$last
  
  # Get mean and sd of ability in sample
  pv = plausible_values_(respData, parms)
  new_mu = mean(pv$PV1)
  new_sigma = sd(pv$PV1)
  
  pgb$tick(n_booklets)
  
  ssI = design %>%
    inner_join(parms$inputs$ssI, by = 'item_id') %>% 
    arrange(.data$booklet_id, .data$first)
  
  a = parms$inputs$ssIS$item_score
  b = parms$est$b
  
  # some common summary statistics
  if (parms$inputs$method == "Bayes")
  {
    n_bayes = nrow(parms$est$b)
    # can do always use all bayes iterations?
    if(n_bayes<nDraws)
    {
      warning('nDraws set to ', n_bayes, '. See the help page for details')
      nDraws = n_bayes
    }
    iter_set = round(seq(1, n_bayes, n_bayes/nDraws))
    iter_set  = iter_set + (n_bayes - iter_set[nDraws])
    ref_theta = sapply(iter_set, function(iter)
    {
      theta_MLE(b[iter,],a, ref_first, ref_last)$theta[pass_fail+1]
    })

  } else
  {
    ref_theta = theta_MLE(b,a,ref_first,ref_last)$theta[pass_fail+1]
  }
  
  bk_results = lapply(split(ssI, ssI$booklet_id), function(dsg)
  {
    scores = 0:sum(a[dsg$last])
    p_new = plausible_scores(respData, parms=parms, items = dsg$item_id) %>%
        count(.data$PS1) %>%
        right_join(tibble(PS1 = as.integer(0:sum(a[dsg$last]))), by='PS1') %>%
        mutate(n = coalesce(.data$n,0L)) %>%
        arrange(.data$PS1) %>%
        pull(.data$n)
    
    p_new = p_new/sum(p_new)
    
    pgb$tick()
    
    max_score = sum(a[dsg$last])
    ref_range = (pass_fail:max_score) + 1
    probs = matrix(0,max_score+1, nDraws)
    
    # some constants for cpp func
    first0 = dsg$first - 1L
    last0 = dsg$last - 1L
    bcni = c(0L,length(first0))
    cbk = integer(length(scores))
    cmu = rep(new_mu,length(scores))
    
    if (parms$inputs$method == "Bayes")
    {
      eq_score = rep(-1L, length(iter_set))
      tel = 1
      # start values
      spv = matrix(theta_jEAP(colMeans(b),a,dsg$first, dsg$last)$theta,ncol=1)
      for (iter in iter_set)
      {
        eq_score[tel] = min(which(theta_MLE(b[iter,],a, dsg$first, dsg$last)$theta >= ref_theta[tel])) - 1
        #spv = pv_recycle_sorted(b[iter,], a,dsg$first, dsg$last, scores, npv=1, mu=new_mu, sigma=new_sigma)
        PV_sve(b[iter,],a, first0, last0, bcni,cbk, scores, cmu, new_sigma,spv, niter=50L)
        
        prf = pscore(spv, b[iter,],a,ref_first,ref_last)
        probs[,tel] = apply(prf, 2, function(x) sum(x[ref_range]))
        tel=tel+1
        pgb$tick()
      }
    } else
    {
      eq_score = min(which(theta_MLE(b,a,dsg$first, dsg$last)$theta >= ref_theta))-1 
      
      #spv = pv_recycle_sorted(b, a, dsg$first, dsg$last, scores, npv=nDraws, mu=new_mu, sigma=new_sigma)
      # start values
      spv = matrix(theta_jEAP(b,a,dsg$first, dsg$last)$theta,ncol=1)
      
      for (iter in 1:nDraws)
      {
        PV_sve(b,a, first0, last0, bcni,cbk, scores, cmu, new_sigma,spv, niter=50L)
        prf = pscore(spv, b, a, ref_first, ref_last)
        #prf = pscore(spv[,iter], b, a, ref_first, ref_last)
        probs[,iter] = apply(prf, 2, function(x) sum(x[ref_range]))
        pgb$tick()
      }
    }
    
    
    p_pass_given_new = rowMeans(probs) 
    
    ## additional
    tp = rev(cumsum(rev(p_pass_given_new*p_new)))/rev(cumsum(rev(p_new)))
    tp_rate = rev(cumsum(rev(p_pass_given_new*p_new)))/tp[1]
    tn = rev(cumsum(rev((1-p_pass_given_new)*p_new)))/rev(cumsum(rev(p_new)))
    fp_rate =  rev(cumsum(rev((1-p_pass_given_new)*p_new)))/tn[1]
    

    list(prob_to_pass = 
           tibble(score_new = 0:max_score, probability_to_pass = p_pass_given_new, true_positive = tp,
                  sensitivity = tp_rate, specificity = 1-fp_rate,
                  true_positive_rate = tp_rate, false_positive_rate = fp_rate, proportion=p_new),
         equated_score=eq_score)
  })

  
  out = list(est = bk_results, 
             inputs = list(method = parms$inputs$method, 
                           ref_items = ref_items, pass_fail = pass_fail,
                           design = design)) 
  class(out) = append('p2pass', class(out))
  out
}






print.p2pass = function(x,...)
{
  cat(paste('Equating information for', length(x$est),
              'booklets. Use `coef()` to extract the statistics or `plot()` to',
              'plot the probability to pass, sensitivity/specificity or a ROC-curve.\n'))
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
  lapply(object$est, function(bk)
    {
      select(bk$prob_to_pass,!ends_with('rate'))
    }
  ) %>%
    bind_rows(.id = 'booklet_id') %>%
    df_format()
}


##########################################
#' A plot method for probability_to_pass
#'
#' Plot equating information from probability_to_pass 
#'
#'
#' @param x An object produced by function \code{\link{probability_to_pass}}
#' @param booklet_id vector of booklet_id's to plot, if NULL all booklets are plotted
#' @param what information to plot, 'equating', 'sens/spec', 'roc, or 'all'
#' @param ... Any additional plotting parameters; e.g., cex = 0.7.
#' @method plot p2pass
#'
plot.p2pass = function(x, what = c('all','equating','sens/spec', 'roc'), booklet_id=NULL, ... )
{
  user.args = list(...)
  what = match.arg(what)
  
  if(!is.null(booklet_id))
  {
    bk_results = x$est[booklet_id] # as list on purpose
  } else
  {
    bk_results = x$est
  }
  
  for(bkl in names(bk_results))
  {
    if(what %in% c('all','equating'))
    {
      dflt = list(xlab = paste("score on", bkl), ylab = "Probability to pass", main = "Equating",
                 pch = 16, bty='l', col='grey40')
      ovr = list(x = bk_results[[bkl]]$prob_to_pass$score_new, 
                 y = bk_results[[bkl]]$prob_to_pass$probability_to_pass)
      args = merge_arglists(user.args, dflt, ovr)
      do.call(plot, args)
      if( x$inputs$method == 'Bayes')
      {
        breaks = seq(min(bk_results[[bkl]]$equated_score)-.5, max(bk_results[[bkl]]$equated_score)+.5)
        hist(bk_results[[bkl]]$equated_score, add=TRUE, freq=FALSE, breaks=breaks)
        do.call(points,args)
      } else
      {
        abline(v=bk_results[[bkl]]$equated_score)
      }
      
    }
    if(what %in% c('all','sens/spec'))
    {
      dflt = list(xlab = paste("score on", bkl),main="Sensitivity/Specificity", ylab = "sens./spec.",
                  pch = 16, ylim=c(0,1), bty='l',col='#1a9641')
      ovr = list(x = bk_results[[bkl]]$prob_to_pass$score_new, 
                 y = bk_results[[bkl]]$prob_to_pass$sensitivity,
                 type="o")
      do.call(plot, merge_arglists(user.args, dflt, ovr))
      
      
      lnargs = dropNulls(merge_arglists(
        user.args[c('lty','lwd','col','pch')],
        list(col='#d7191c',pch = 16),
        list(type="o",
             x = bk_results[[bkl]]$prob_to_pass$score_new,
             y = bk_results[[bkl]]$prob_to_pass$specificity)))
        
      do.call(lines, lnargs)
    }
    if(what %in% c('all','roc'))
    {
      dflt = list(xlab = '1 - specificity', ylab='sensitivity', main="ROC", 
                  xlim=c(0,1), ylim=c(0,1), bty='l',col='#1a9641')
      ovr = list(x = 1-bk_results[[bkl]]$prob_to_pass$specificity, 
                 y = bk_results[[bkl]]$prob_to_pass$sensitivity,
                 type="l")
      do.call(plot, merge_arglists(user.args, dflt, ovr))
      text(ovr$x, ovr$y, bk_results[[bkl]]$prob_to_pass$score_new, cex=0.7, offset = 0)
      abline(0,1,lty=2, col="grey")
    }
  }
  
  invisible(NULL)
}

