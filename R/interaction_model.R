
#' Estimate the Interaction and the Rasch model
#'
#' Estimate the parameters of the Interaction model and the Rasch model
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @return An object of class \code{rim} holding results
#' for the Rasch model and the interaction model.
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function fits the
#' Rasch model and the interaction model on a complete rectangular array of
#' responses. This typically consist of responses to items in one booklet but can also consist of
#' the intersection (common items) of two or more booklets. If the intersection is empty
#' (no common items for all persons), the function will exit with an error message.
#'
#' @seealso \code{\link{plot.rim}}, \code{\link{fit_domains}}
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#'
#' m = fit_inter(db, booklet_id=='agg')
#' plot(m, "S1DoScold", show.observed=TRUE)
#'
#' close_project(db)
#' }
#'
#'
fit_inter <- function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))

  respData = get_resp_data(dataSrc, qtpredicate, env = caller_env())

  if(nrow(respData$x)==0) stop('no responses to analyse')

  # make sure we have an intersection
  if(length(unique(respData$design$booklet_id)) > 1)
  {
    items = Reduce(intersect, split(respData$design$item_id, respData$design$booklet_id))
    if(length(items)==0) stop(paste('The intersection of responses in your data is empty.',
                                    'The interaction model cannot be computed concurrently for a whole design of test forms.',
                                    'See the help for more information'))

    respData$design = tibble(booklet_id='b', item_id = items)

    respData$x = respData$x %>%
      semi_join(respData$design, by='item_id') %>%
      group_by(.data$person_id) %>%
      mutate(sumScore = sum(.data$item_score), booklet_id = 'b') %>%
      ungroup()
  }

  if(is.null(qtpredicate)){ grpName = 'All'} else{ grpName = as.character(qtpredicate)}

  # add 0 scores if necessary
  ssScoreLev = respData$x %>%
    group_by(.data$item_id, .data$item_score) %>%
    summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>%
    ungroup() %>%
    full_join(tibble(item_id=respData$design$item_id, item_score=0), by = c("item_id","item_score")) %>%
    arrange(.data$item_id, .data$item_score)

  ssScoreLev[is.na(ssScoreLev)] = 0

  ssItemLev = ssScoreLev %>%
    group_by(.data$item_id) %>%
    summarise(nCat = n(), N = sum(.data$sufI), sufC = sum(.data$sufC)) %>%
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1, last = cumsum(.data$nCat))  %>%
    ungroup() %>%
    arrange(.data$item_id)

  plt = respData$x %>%
    group_by(.data$item_id, .data$sumScore) %>%
    summarise(meanScore = mean(.data$item_score), N = n()) %>%
    ungroup()

  maxTotScore = sum(tapply(ssScoreLev$item_score, ssScoreLev$item_id, max))

  totalLev = plt[plt$item_id==plt$item_id[1],] %>%
    right_join(tibble(sumScore=0:maxTotScore), by="sumScore")

  totalLev$N[is.na(totalLev$N)] = 0

  ss = list(group=grpName, il=ssItemLev, sl=ssScoreLev, tl=totalLev, plt=plt)

  # changed because I don't like this, if estimation throws an error,
  # this function should throw one too, not return a mostly useless object
  #result = try(EstIM(ss))
  #if (inherits(result, "try-error")) result=NULL

  result = EstIM(ss)

  # add the regressions, convenient for plotting
  #if (!is.null(result))
  #{
    C = rep(1:nrow(ss$il), ss$il$nCat)
    ctrRM=ittotmat(result$bRM, result$cRM[C], ss$sl$item_score, ss$il$first, ss$il$last)
    ctrIM=ittotmat(result$bIM, result$cIM[C], ss$sl$item_score, ss$il$first, ss$il$last)
    mm = sweep(model.matrix(~0+ss$sl$item_id), 1, ss$sl$item_score, '*')
    itrRM = as.data.frame(crossprod(mm, ctrRM))
    itrIM = as.data.frame(crossprod(mm, ctrIM))
    row.names(itrRM) = row.names(itrIM) = ss$il$item_id
    regs = list(ctrRM, ctrIM, itrRM, itrIM)
    names(regs) = c('ctrRM','ctrIM','itrRM','itrIM')
  #} else {
  #  regs = NULL
  #}
  outpt = list(est=result, ss=ss, regs=regs)
  class(outpt) = append("rim",class(outpt))
  outpt
}






#' Estimate the Rasch and the Interaction model per domain
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
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
#' @seealso \code{\link{plot.rim}}, \code{\link{fit_inter}}, \code{\link{add_item_properties}}
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db")
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' mSit = fit_domains(db, item_property= "situation")
#' plot(mSit)
#'
#' close_project(db)
#' }
#'
fit_domains = function(dataSrc, item_property, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))

  respData = get_resp_data(dataSrc, qtpredicate, extra_columns=item_property, env=caller_env())
  if(nrow(respData$x) == 0) stop('no data to analyse')

  if(length(unique(respData$design$booklet_id)) > 1)
  {
    items = Reduce(intersect, split(respData$design$item_id, respData$design$booklet_id))
    if(length(items)==0) stop(paste('The intersection of responses in your data is empty.',
                                    'The interaction model cannot be computed concurrently for a whole design of test forms.',
                                    'See the help for more information'))

    respData$design = tibble(booklet_id='b', item_id = items)

    respData$x = respData$x %>%
      semi_join(respData$design, by='item_id')
  }

  respData$x = respData$x %>%
    group_by(.data$person_id, .data[[!!item_property]]) %>%
    summarise(item_score=sum(.data$item_score)) %>%
    ungroup() %>%
    group_by(.data$person_id) %>%
    mutate(sumScore = sum(.data$item_score)) %>%
    ungroup() %>%
    rename(item_id = .data[[!!item_property]]) %>%
    add_column(booklet_id='b')

  respData$design = tibble(booklet_id = 'b', item_id = unique(respData$x$item_id))

  fit_inter(respData)
}

draw_curtains = function(qnt)
{
 if(!is.null(qnt))
 {
   usr = graphics::par('usr')
   graphics::rect(usr[1], usr[3], qnt[1], usr[2], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
   graphics::rect(qnt[2], usr[3], usr[2], usr[4], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
   # rect will sometimes cover the axis (redrawing the axis would probably solve this too)
   abline(h=usr[3])
   abline(v=usr[1])
 }
}

##########################################
#' A plot method for the interaction model
#'
#' Plot the item-total regressions fit by the interaction (or Rasch) model
#'
#'
#' @param x An object produced by function \code{fit_inter}
#' @param items The items to plot (item_id's). If NULL, all items will be plotted
#' @param summate If FALSE, regressions for polytomous items will be shown for each
#' response option separately; default is TRUE.
#' @param overlay If TRUE and more than one item is specified, there will be two plots,
#' one for the Rasch model and the other for the interaction model, with all items
#' overlayed; otherwise, one plot for each item with the two models overlayed. Ignored
#' if summate is FALSE. Default is FALSE
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param show.observed If TRUE, the observed proportion correct at each sum score
#' will be shown as dots. Default is FALSE.
#' @param ... Any additional plotting parameters.
#' @details
#' Customisation of title and subtitle can be done by using the arguments main and sub.
#' These arguments can contain references to the variables item_id (if overlay=FALSE) or model (if overlay=TRUE)
#' by prefixing them with a dollar sign, e.g. plot(m, main='item: $item_id')
#' @method plot rim
#'
plot.rim <- function(x, items=NULL, summate=TRUE, overlay=FALSE,
                     nc=1, nr=1, curtains=10, show.observed=FALSE, ...){
  allItems = x$ss$il$item_id
  if(!is.null(items))
  {
    if(length(setdiff(items, allItems)) > 0) stop(paste('item(s):', paste0(setdiff(items, allItems), collapse=', '), 'not found'))
  } else
  {
    items = allItems
  }
  user.args = list(...)

  qua = curtains/200
  qnt=NULL
  if(qua>0 && qua<.5) {
    qnt = quantile(rep(as.integer(x$ss$tl$sumScore), x$ss$tl$N), c(qua,1-qua))
  }

  npic = length(items)
  if (length(items)==1) nr=nc=1
  if (overlay & !summate) overlay=FALSE
  if (overlay) {
    # only summate possible
    if (nr*nc==2) graphics::layout(matrix(1:2,nr,nc)) else graphics::layout(1)
    # do the Rasch model
    #
    z = x$regs$itrRM
    z = z[row.names(z) %in% items,]
    maxScore = ncol(z)-1

    plot.args = merge_arglists(user.args,
                               default=list(main="$model model",xlab="Test score",
                                            ylab="Item score", bty='l'),
                               override=list(x = c(0,maxScore),y= c(0,max(z)), type="n"))

    plot.args$main = fstr(plot.args$main, list(model='Rasch'))
    plot.args$sub = fstr(plot.args$sub, list(model='Rasch'))

    do.call(graphics::plot, plot.args)
    draw_curtains(qnt)
    
    for (i in 1:npic) graphics::lines(0:maxScore,z[i,]) # the actual lines
    lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      graphics::points(lx[i], z[i,lx[i]+1], co="white", cex=1.6, pch=19)
      graphics::text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    # do the Interaction model
    #
    z = x$regs$itrIM
    z = z[row.names(z) %in% items,]
    maxScore = ncol(z)-1

    plot.args = merge_arglists(user.args,
                               default=list(main="$model model",xlab="Test score",
                                            ylab="Item score", bty='l'),
                               override=list(x = c(0,maxScore),y= c(0,max(z,na.rm=TRUE)), type="n"))

    plot.args$main = fstr(plot.args$main, list(model='Interaction'))
    plot.args$sub = fstr(plot.args$sub, list(model='Interaction'))

    do.call(graphics::plot, plot.args)
    draw_curtains(qnt)
    
    for (i in 1:npic) graphics::lines(0:maxScore,z[i,]) # the actual lines
    lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      graphics::points(lx[i], z[i,lx[i]+1], co="white", cex=1.6, pch=19)
      graphics::text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    
    # end of overlay

  } else {
    # not overlay: do many plots
    ly = my_layout(npic, nr, nc)
    graphics::layout(matrix(1:(ly$nr*ly$nc), byrow=TRUE, ncol=ly$nc))
    if (summate) {
      # for each item in turn, do both models (with summation), and plot
      zI = x$regs$itrIM
      zR = x$regs$itrRM
      maxScore = ncol(zR)-1
      for (i in items) {
        mxY = max(zR[row.names(zR)==i,],na.rm=TRUE)

        plot.args = merge_arglists(user.args,
                                   default=list(main='Item $item_id',xlab="Test score",
                                                ylab="Item score", bty='l'),
                                   override=list(x = c(0,maxScore),y = c(0,mxY), type="n"))

        plot.args$main = fstr(plot.args$main, list(item_id=i))
        plot.args$sub = fstr(plot.args$sub, list(item_id=i))

        do.call(graphics::plot, plot.args)

        draw_curtains(qnt)
        if(show.observed) {
          plt = x$ss$plt[x$ss$plt$item_id==i,]
          graphics::points(plt$sumScore,plt$meanScore,col="coral",pch=20)
        }
        graphics::lines(0:maxScore, zI[row.names(zI)==i,], col="gray80", lwd=3)
        graphics::lines(0:maxScore, zR[row.names(zR)==i,])
      }
      
    } else {
      zI = x$regs$ctrIM
      zR = x$regs$ctrRM
      maxScore = ncol(zR)-1
      # for each item in turn, similar but possibly multiline and coloured
      for (i in items) {
        prb = zI[x$ss$il$first[x$ss$il$item_id==i]:x$ss$il$last[x$ss$il$item_id==i],]
        pte = bty(nrow(prb), alpha=.6)
        #pte = qcolors(nrow(prb))

        plot.args = merge_arglists(user.args,
                                   default=list(main='Item $item_id',xlab="Test score",
                                                ylab="Probability",bty='l'),
                                   override=list(x = c(0,maxScore),y = 0:1, type="n"))

        plot.args$main = fstr(plot.args$main, list(item_id=i))
        plot.args$sub = fstr(plot.args$sub, list(item_id=i))

        do.call(graphics::plot, plot.args)

        draw_curtains(qnt)
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j], lwd=3)
        }
        prb = zR[x$ss$il$first[x$ss$il$item_id==i]:x$ss$il$last[x$ss$il$item_id==i],]
        pte = bty(nrow(prb))
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j])
        }
      } # eol items
      
    } # eo not summate
    if(nc>1 || nr>1) par(mfrow=c(1,1))
  } # eo not overlay
}



###########################################
#' A print method for the interaction model
#'
#' Print the available items for plots of the Rasch and the interaction models
#'
#'
#' @param x An object produced by function \code{fit_inter}
#' @param ... Further arguments to print
#' @method print rim
#'
print.rim <- function(x, ...){
  do.call(print.default, merge_arglists(list(...), override = list(x=pull(x$ss$il,"item_id"))))
}


###########################################
#' Coefficients for the interaction model
#'
#' Extract parameters for the interaction and the Rasch model
#'
#'
#' @param object An object produced by function \code{fit_inter}
#' @param ... Further arguments currently ignored
#' @method coef rim
#'
coef.rim = function(object, ...) 
{
  #can also add an argument like what=c('parameters','fit')
  x = object
  first=x$ss$il$first
  dRM=x$est$bRM[-first]
  dRM=-log(dRM/dRM[1])
  dIM=x$est$bIM[-first]
  dIM=-log(dIM/dIM[1])
  
  IS = tibble(item_id = x$ss$sl$item_id[-first], item_score = x$ss$sl$item_score[-first],
              delta_rasch = dRM, delta_IM = dIM)
  I = tibble(item_id = x$ss$il$item_id, sigma = log(x$est$cIM), SE_sigma= x$est$se.c, fit_IM=x$est$fit.stats)
  
  inner_join(IS,I,by='item_id') %>% arrange(.data$item_id, .data$item_score)
}

