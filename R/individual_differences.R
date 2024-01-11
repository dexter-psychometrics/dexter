##########################################
#' Test individual differences
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data are used.
#' @details This function uses a score distribution to test whether there are individual 
#' differences in ability. First, it estimates ability based on the score distribution. Then, 
#' the observed distribution is compared to the one expected from the single estimated ability.
#' The data are typically from one booklet but can also consist of 
#' the intersection (i.e., the common items) of two or more booklets. If the intersection is empty
#' (i.e., no common items for all persons), the function will exit with an error message.
#' @return An object of type tind. Printing the object will show test results. 
#' Plotting it will produce a plot of expected and observed score frequencies. 
#' The former under the hypothesis that there are no individual differences.
#'
#'@examples
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules, ":memory:")
#' add_booklet(db, verbAggrData, "agg")
#' 
#' dd = individual_differences(db)
#' print(dd)
#' plot(dd)
#' 
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
individual_differences = function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  check_dataSrc(dataSrc)
  
  respData = get_resp_data(dataSrc, qtpredicate, env = env) |>
	  intersection_rd()
  
  parms = fit_enorm(respData)
  b = parms$est$b
  a = parms$inputs$ssIS$item_score
  
  # intersection, so just one booklet
  first = parms$inputs$ssI$first
  last = parms$inputs$ssI$last
  observed = parms$inputs$scoretab$N
  m = sum(observed)
  lambda = parms$est$lambda$lambda
  observed_smooth = suppressWarnings({ENORM2ScoreDist(b,a,lambda,first,last)$n.smooth})

  theta.est = theta_score_distribution(b,a,first,last,observed)
  expected = pscore(theta.est,b,a,first,last)[,1,drop=TRUE]
  chi = chisq.test(x=observed,p=expected,simulate.p.value = TRUE)
  
  
  inputs = list(items=parms$inputs$ssI$item_id, m = m, max.score=sum(a[last]), 
                observed=observed, observed_smooth=observed_smooth, xpr=as.character(qtpredicate))
  est =list(test=chi, theta=theta.est)
  outpt = list(inputs=inputs, est=est)
  class(outpt) = append("tind",class(outpt))
  outpt
}



print.tind=function(x,...)
{
  cat("Chi-Square Test for the hypothesis that respondents differ in ability:\n")
  print(x$est$test,...)
}


coef.tind=function(object,...)
{
	object$est
}

plot.tind=function(x,...)
{
  user.args = list(...); leg.args = list()
  col = c("#4DAF4A","gray")
  if('col' %in% names(user.args))
  {
    col = coalesce(as.character(user.args$col[1:2]),col)
    user.args$col = NULL
  }
  
  if(length(names(user.args))>0)
  {
    leg.args = user.args[endsWith(names(user.args),'legend')]
    names(leg.args) = gsub('\\.legend$','',names(leg.args))
    user.args = user.args[!endsWith(names(user.args),'legend')]
  }
  
  mx.frq=max(max(x$est$test$expected),max(x$est$test$observed))
  default.args = list(ylim=c(0,mx.frq),
                      xlab="Test-score", 
                      ylab="Frequency",
                      cex=0.7, 
                      bty='l',
                      pch=19)
  override.args = list(x=0:x$inputs$max.score, y=x$est$test$observed,col=col[1])
 
  do.call(plot,merge_arglists(user.args,override=override.args, default=default.args))
  
  lines(0:x$inputs$max.score,x$inputs$observed_smooth,col=col[1])
  lines(0:x$inputs$max.score,x$est$test$expected,col=col[2],pch=19,cex=0.7)
  
  do.call(graphics::legend,
            merge_arglists(leg.args, 
                           default=list(x="topleft", bty = "n",
                                        lwd = 1, cex = 0.7, col = col, pch = c(19,NA),inset=0.01,
                                        legend=c("observed", "expected")),
                           override=list(lty = 1, col=col)))

  
  invisible(NULL)
}