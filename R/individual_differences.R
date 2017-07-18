##########################################
#' Test individual differences
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @details This function uses a score distribution to test whether there are individual 
#' differences in ability. First, it estimates ability based on the score distribution. Then, 
#' the observed distribution is compared to the one expected from the single estimated ability.
#' The data are typically from one booklet but can also consist of 
#' the intersection (i.e., the common items) of two or more booklets. If the intersection is empty 
#' (no common items for all persons), the function will exit with an error message.
#' @return Print will show test results. Plot will produce a plot of expected and observed score frequencies.
#'
#'@examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' dd = individual_differences(db)
#' print(dd)
#' plot(dd)
#' }
#' 
individual_differences <- function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  respData = get_responses_(dataSrc,qtpredicate, env=caller_env()) 
  
  if (nrow(respData)<1) stop("No responses to analyse")
  respData = 
    respData %>%
    spread_(key_col='item_id', value_col='item_score') %>%
    select(which(c(TRUE, !is.na(colSums(.[,-1]))))) 
  
  if(ncol(respData)==1) stop (paste('No common items across booklets.'))
  
  itcol = setdiff((names(respData)), 'person_id')
  respData = respData %>%
    gather_(key_col='person_id', value_col='item_score', itcol)
  colnames(respData)[2]="item_id"  ### this is ugly
  
  parms=fit_enorm(respData)
  b=parms$est$b
  a=parms$inputs$ssIS$item_score
  first=parms$inputs$bkList[[1]]$first
  last=parms$inputs$bkList[[1]]$last
  m=parms$inputs$bkList[[1]]$m
  
  observed=parms$inputs$bkList[[1]]$scoretab
  theta.est = theta_score_distribution(b,a,first,last,observed)
  expected = pscore(theta.est,b,a,first,last)
  max.score=sum(a[last])
  
  chi = chisq.test(x=observed,p=expected,simulate.p.value = TRUE)
  
  inputs = list(items=parms$inputs$bkList[[1]]$items, m=m, max.score=max.score, observed=observed, xpr=as.character(qtpredicate))
  est =list(test=chi, theta=theta.est)
  outpt = list(inputs=inputs, est=est)
  class(outpt) = append("tind",class(outpt))
  outpt
}

print.tind=function(x,...)
{
  print("Chi-Square Test for the hypothesis that all respondents have the same ability:")
  print(x$est$test,...)
}

plot.tind=function(x,...)
{
  user.args = list(...)
  mx.frq=max(max(x$est$test$expected),max(x$est$test$observed))
  default.args = list(col="green",
                      ylim=c(0,mx.frq),
                      xlab="Test-score", 
                      ylab="Frequency",
                      cex=0.7, 
                      pch=19)
  override.args = list(x=0:x$inputs$max.score, y=x$est$test$observed)
  
  do.call(plot,merge_arglists(user.args,override=override.args, default=default.args))

  lines(0:x$inputs$max.score,x$est$test$expected,col="gray",pch=19,cex=0.7)
  legend("topleft", legend = c("observed", "expected"), bty = "n",
         lwd = 1, cex = 0.7, col = c("green", "gray"), lty = c(NA,1), pch = c(19,NA))
}