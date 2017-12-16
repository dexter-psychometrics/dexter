##########################################
#' Test individual differences
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data are used.
#' @param degree The degree of a polynomial used to smooth observed score distribution. 
#' @details This function uses a score distribution to test whether there are individual 
#' differences in ability. First, it estimates ability based on the score distribution. Then, 
#' the observed distribution is compared to the one expected from the single estimated ability.
#' The data are typically from one booklet but can also consist of 
#' the intersection (i.e., the common items) of two or more booklets. If the intersection is empty 
#' (no common items for all persons), the function will exit with an error message.
#' @return an object of type tind. Printing the object  will show test results. 
#' Plotting it will produce a plot of expected and observed score frequencies. 
#' The former under the hypothesis that there are no individual differences.
#'
#'@examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' dd = individual_differences(db)
#' print(dd)
#' plot(dd)
#' 
#' close_project(db)
#' }
#' 
individual_differences <- function(dataSrc, predicate = NULL, degree=7)
{
  qtpredicate = eval(substitute(quote(predicate)))
  
  respData = get_resp_data(dataSrc, qtpredicate, env = caller_env()) 
  
  # make sure we have an intersection
  if(length(unique(respData$design$booklet_id)) > 1)
  {
    respData$design = tibble(booklet_id='b', 
                             item_id = Reduce(intersect, split(respData$design$item_id, respData$design$booklet_id)))
    
    respData$x = respData$x %>%
      semi_join(respData$design, by='item_id') %>%
      group_by(.data$person_id) %>%
      mutate(sumScore = sum(.data$item_score), booklet_id = 'b') %>%
      ungroup()
  }

  if (nrow(respData$x)<1) stop ('No common items across booklets.')
  
  parms=fit_enorm(respData)
  b=parms$est$b
  a=parms$inputs$ssIS$item_score
  first=parms$inputs$bkList[[1]]$first
  last=parms$inputs$bkList[[1]]$last
  m=parms$inputs$bkList[[1]]$m
  observed=parms$inputs$bkList[[1]]$scoretab
  observed_smooth=ENORM2ScoreDist(parms, degree,booklet_id=parms$inputs$bkList[[1]]$booklet)$n.smooth
  
  theta.est = theta_score_distribution(b,a,first,last,observed)
  expected = pscore(theta.est,b,a,first,last)
  chi = chisq.test(x=observed,p=expected,simulate.p.value = TRUE)
  
  
  inputs = list(items=parms$inputs$bkList[[1]]$items, m=m, max.score=sum(a[last]), 
                observed=observed, observed_smooth=observed_smooth, xpr=as.character(qtpredicate))
  est =list(test=chi, theta=theta.est)
  outpt = list(inputs=inputs, est=est)
  class(outpt) = append("tind",class(outpt))
  outpt
}


print.tind=function(x,...)
{
  cat("Chi-Square Test for the hypothesis that all respondents have the same ability:\n")
  print(x$est$test,...)
}


plot.tind=function(x,...)
{
  user.args = list(...)
  mx.frq=max(max(x$est$test$expected),max(x$est$test$observed))
  default.args = list(col="#4DAF4A",
                      ylim=c(0,mx.frq),
                      xlab="Test-score", 
                      ylab="Frequency",
                      cex=0.7, 
                      bty='l',
                      pch=19)
  override.args = list(x=0:x$inputs$max.score, y=x$est$test$observed)
 
  do.call(plot,merge_arglists(user.args,override=override.args, default=default.args))

  lines(0:x$inputs$max.score,x$est$test$expected,col="gray",pch=19,cex=0.7)
  lines(0:x$inputs$max.score,x$inputs$observed_smooth,col="lightgreen")
  legend("topleft", legend = c("observed", "expected"), bty = "n",
         lwd = 1, cex = 0.7, col = c("#4DAF4A", "gray"), lty = c(NA,1), pch = c(19,NA))
}