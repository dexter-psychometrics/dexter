

# to~do: example    
  

#' The probability to pass on a reference test given a score on a new booklet
#' 
#' Given response data that form a connected design,
#' compute the probability to pass on the reference set conditional on each score on one or more target tests.
#' 
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and  beta. 
#' If uncertainty about parameter estimation should be included
#' in the computations, use a `parms` object computed with `method='Bayes'` and nDraws equal or larger than nDraws in probability_to_pass
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
#' 
probability_to_pass = function(dataSrc, parms, ref_items, pass_fail, predicate = NULL, 
                               target_booklets = NULL, nDraws = 1000)
{
  check_dataSrc(dataSrc)
  check_num(pass_fail, 'integer', .min=0, .length=1)
  check_df(target_booklets, 'item_id', nullable=TRUE)
  
  if(inherits(ref_items,'data.frame') && 'item_id' %in% colnames(ref_items))
  {
    if('booklet_id' %in% colnames(ref_items) && n_distinct(ref_items$booklet_id)>1)
      stop_('ref_itemsshould contain a single booklet')
    ref_items = ref_items$item_id
  }
  if(is.factor(ref_items)) ref_items = as.character(ref_items)
  
  check_character(ref_items)
  
  
  
  pb = get_prog_bar(nsteps=nDraws+10,retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env)
  
  if(is.null(target_booklets))
  {
    design = respData$design
  } else
  {
    design = target_booklets %>%
      mutate(item_id = factor(as.character(.data$item_id), levels = levels(respData$design$item_id)))
        
    if(anyNA(design$item_id))
      stop('One or more items in design does not occur in your dataSrc')
    
    if(!('booklet_id' %in% colnames(design))) 
      design$booklet_id = 'target items'
  }
  
  
  target_parms = simplify_parms(parms, design=design, draw='sample') 
  
  ref_items = tibble(item_id = ffactor(unique(ref_items), levels=levels(respData$design$item_id)))
  
  if(anyNA(ref_items$item_id))
    stop('One or more of your reference items does not occur in your dataSrc')
  
  ref_parms = simplify_parms(parms,design=ref_items)
 
  n_booklets = n_distinct(design$booklet_id)
  pb$set_nsteps(10+(nDraws+10)*n_booklets)
  

  # Get mean and sd of ability in sample
  pb$new_area(10)
  pv = plausible_values_(respData, parms, parms_draw='average')$pv
  new_mu = mean(pv$PV1)
  new_sigma = sd(pv$PV1)
  
  max_ref_score = sum(ref_parms$a[ref_parms$design$last])
  ref_range = (pass_fail:max_ref_score) + 1L
  pb$close_area()
  
  # some common summary statistics for the reference items
  if (ref_parms$method == "Bayes")
  {
    n_bayes = nrow(ref_parms$b)
    # can do always use all bayes iterations?
    if(n_bayes < nDraws)
    {
      warning('nDraws set to ', n_bayes, '. See the help page for details')
      nDraws = n_bayes
    }
    iter_set = round(seq(1, n_bayes, n_bayes/nDraws))
    iter_set  = iter_set + (n_bayes - iter_set[nDraws])
    ref_theta = sapply(iter_set, function(iter)
    {
      theta_MLE(ref_parms$b[iter,],ref_parms$a, ref_parms$design$first, ref_parms$design$last)$theta[pass_fail+1]
    })
    
  } else
  {
    ref_theta = theta_MLE(ref_parms$b,ref_parms$a, ref_parms$design$first, ref_parms$design$last)$theta[pass_fail+1]
  }
  
  max_cores = get_ncores(desired = 32L, maintain_free = 1L)
  
  bk_results = lapply(split(target_parms$design, target_parms$design$booklet_id), function(dsg)
  {
    pb$new_area(10)
    max_score = sum(target_parms$a[dsg$last])
    scores = 0:max_score
    ps_new = plausible_scores(respData, parms=parms, items = dsg$item_id) %>%
      count(.data$PS1) %>%
      right_join(tibble(PS1 = scores), by='PS1') %>%
      mutate(n = coalesce(.data$n,0L)) %>%
      arrange(.data$PS1) %>%
      pull(.data$n)
    
    pb$close_area()
    
    ps_new = ps_new/sum(ps_new)
    probs = matrix(0,max_score+1, nDraws)
    
    # some constants for cpp func
    bcni = c(0L,length(dsg$first0))
    cbk = integer(length(scores))
    cmu = rep(new_mu,length(scores))
    
    
    if(target_parms$method == "Bayes")
    {
      eq_score = rep(-1L, length(iter_set))
      tel = 1
      # start values
      spv = matrix(theta_EAP_GH(colMeans(target_parms$b),target_parms$a,dsg$first, dsg$last, mu=new_mu, sigma=new_sigma)$theta,ncol=1)
      for(iter in iter_set)
      {
        eq_score[tel] = min(which(theta_MLE(target_parms$b[iter,],target_parms$a, dsg$first, dsg$last)$theta >= ref_theta[tel])) - 1

        PV_sve(target_parms$b[iter,],target_parms$a, dsg$first0, dsg$last0, bcni,cbk, scores, cmu, new_sigma,spv, niter=30L, max_cores=max_cores)
        
        prf = pscore(spv, ref_parms$b[iter,],ref_parms$a,ref_parms$design$first,ref_parms$design$last)
        probs[,tel] = apply(prf, 2, function(x) sum(x[ref_range]))
        tel=tel+1
        pb$tick()
      }
    } else
    {
      eq_score = min(which(theta_MLE(target_parms$b,target_parms$a,dsg$first, dsg$last)$theta >= ref_theta))-1 
      
      # start values
      spv = matrix(theta_EAP_GH(target_parms$b,target_parms$a,dsg$first, dsg$last, mu=new_mu, sigma=new_sigma)$theta,ncol=1)
      
      for (iter in 1:nDraws)
      {
        PV_sve(target_parms$b,target_parms$a, dsg$first0, dsg$last0, bcni,cbk, scores, cmu, new_sigma,spv, niter=30L, max_cores=max_cores)
        prf = pscore(spv, ref_parms$b, ref_parms$a, ref_parms$design$first, ref_parms$design$last)
        probs[,iter] = apply(prf, 2, function(x) sum(x[ref_range]))
        pb$tick()
      }
    }

    p_pass_given_new = rowMeans(probs) 
    
    ## additional
    tp = rev(cumsum(rev(p_pass_given_new*ps_new)))/rev(cumsum(rev(ps_new)))
    tp_rate = rev(cumsum(rev(p_pass_given_new*ps_new)))/tp[1]
    tn = rev(cumsum(rev((1-p_pass_given_new)*ps_new)))/rev(cumsum(rev(ps_new)))
    fp_rate =  rev(cumsum(rev((1-p_pass_given_new)*ps_new)))/tn[1]
    
    
    list(prob_to_pass = 
           tibble(score_new = 0:max_score, probability_to_pass = p_pass_given_new, true_positive = coalesce(tp,1),
                  sensitivity = tp_rate, specificity = 1-fp_rate,
                  true_positive_rate = tp_rate, false_positive_rate = fp_rate, proportion=ps_new),
         equated_score=eq_score)
  })
  
  
  out = list(est = bk_results, 
             inputs = list(method = ref_parms$method, 
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
#' \item{true_positive}{proportion that correctly passes}
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

