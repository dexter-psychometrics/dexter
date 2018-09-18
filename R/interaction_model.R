
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
fit_inter = function(dataSrc, predicate = NULL)
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, env = env)

  if(nrow(respData$x)==0) stop('no responses to analyse')

  # test if we have complete data, otherwise make a rectangular matrix by throwing away items
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

  # statistics per item-score combination
  # if the score 0 does not occur for an item, it is added with sufI=0 and sufC=0
  ssIS = respData$x %>%
    group_by(.data$item_id, .data$item_score) %>%
    summarise(sufI=n(), sufC_ = sum(.data$item_score * .data$sumScore)) %>%
    ungroup() %>%
    full_join(tibble(item_id=respData$design$item_id, item_score=0L), by = c("item_id","item_score")) %>%
    mutate(sufI = coalesce(.data$sufI, 0L), sufC_ = coalesce(.data$sufC_,0L)) %>%
    arrange(.data$item_id, .data$item_score)

  # statistics per item
  ssI = ssIS %>%
    group_by(.data$item_id) %>%
    summarise(nCat = n(), sufC = sum(.data$sufC_), item_maxscore = max(.data$item_score)) %>%
    ungroup() %>%
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1L, last = cumsum(.data$nCat))  %>%
    arrange(.data$item_id)

  # mean item score and weight (as n responses) for each individual test score
  # used for plotting
  plt = respData$x %>%
    group_by(.data$item_id, .data$sumScore) %>%
    summarise(meanScore = mean(.data$item_score), N = n()) %>%
    ungroup()

  # theoretical max score on the test
  maxTestScore = sum(ssI$item_maxscore)

  # scoretab from 0:maxscore (including unachieved and impossible scores)
  scoretab = respData$x %>%
    distinct(.data$person_id, .keep_all=TRUE) %>%
    group_by(.data$sumScore) %>%
    summarise(N=n()) %>%
    ungroup() %>%
    right_join(tibble(sumScore=0L:maxTestScore), by="sumScore") %>%
    mutate(N=coalesce(.data$N, 0L)) %>%
    arrange(.data$sumScore)
  

  result = EstIM(first = ssI$first, last = ssI$last, nCat = ssI$nCat, a = ssIS$item_score, 
                 sufI = ssIS$sufI, sufC = ssI$sufC, scoretab = scoretab$N)

  # add the regressions, convenient for plotting

  C = rep(1:nrow(ssI), ssI$nCat)
  
  ctrRM = ittotmat(result$bRM, result$cRM[C], ssIS$item_score, ssI$first, ssI$last)
  ctrIM = ittotmat(result$bIM, result$cIM[C], ssIS$item_score, ssI$first, ssI$last)
  mm = sweep(model.matrix(~0 + ssIS$item_id), 1, ssIS$item_score, '*')
  itrRM = as.data.frame(crossprod(mm, ctrRM))
  itrIM = as.data.frame(crossprod(mm, ctrIM))
  row.names(itrRM) = row.names(itrIM) = ssI$item_id

  output = list(est = result, 
                inputs = list(ssI = ssI, ssIS = ssIS, plt = plt, scoretab=scoretab), 
                regs = list(ctrRM = ctrRM, ctrIM = ctrIM, itrRM = itrRM, itrIM = itrIM))
  class(output) = append("rim", class(output))
  output
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
  env = caller_env()
  respData = get_resp_data(dataSrc, qtpredicate, extra_columns=item_property, env = env)
  if(nrow(respData$x) == 0) stop('no data to analyse')

  # test if we have complete data, otherwise make a rectangular matrix by throwing away items
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

  # adapt the respdata object by making new polytomous items based on the domains
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

  # call fit_inter on the adapted data
  fit_inter(respData)
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
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param show.observed If TRUE, the observed proportion correct at each sum score
#' will be shown as dots. Default is FALSE.
#' @param ... Any additional plotting parameters.
#' @details
#' Customization of title and subtitle can be done by using the arguments main and sub.
#' These arguments can contain references to the variables item_id (if overlay=FALSE) or model (if overlay=TRUE)
#' by prefixing them with a dollar sign, e.g. plot(m, main='item: $item_id')
#' @method plot rim
#'
plot.rim = function(x, items=NULL, summate=TRUE, overlay=FALSE,
                     curtains=10, show.observed=FALSE, ...){
  allItems = x$inputs$ssI$item_id
  if(!is.null(items))
  {
    if(length(setdiff(items, allItems)) > 0) stop(paste('item(s):', paste0(setdiff(items, allItems), collapse=', '), 'not found'))
  } else
  {
    items = allItems
  }
  user.args = list(...)

  if(any(c('nr','nc') %in% names(user.args)))
  {
    warning('Arguments `nr` and `nc` have been deprecated and are ignored, you can achieve the desired ',
            'result by using par(mfrow=c(nr,nc)).')
    user.args[['nc']] = NULL
    user.args[['nr']] = NULL
  }
  
  qua = curtains/200
  qnt=NULL
  if(qua>0 && qua<.5) {
    qnt = quantile(rep(as.integer(x$inputs$scoretab$sumScore), x$inputs$scoretab$N), c(qua,1-qua))
  }

  npic = length(items)

  if (overlay && !summate) overlay=FALSE
  if (overlay) {
    # only summate possible
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

    do.call(plot, plot.args)
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

    do.call(plot, plot.args)
    draw_curtains(qnt)
    
    for (i in 1:npic) graphics::lines(0:maxScore, z[i,]) # the actual lines
    lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      graphics::points(lx[i], z[i,lx[i]+1], co="white", cex=1.6, pch=19)
      graphics::text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    
    # end of overlay

  } else {
    # not overlay: do many plots
    
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

        do.call(plot, plot.args)

        draw_curtains(qnt)
        if(show.observed) {
          plt = filter(x$inputs$plt, .data$item_id == i)
          graphics::points(plt$sumScore, plt$meanScore, col="coral",pch=20)
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
        ssI = filter(x$inputs$ssI, .data$item_id==i)
        prb = zI[ssI$first : ssI$last, ]
        
        pte = bty(nrow(prb), alpha=.6)
        #pte = qcolors(nrow(prb))

        plot.args = merge_arglists(user.args,
                                   default=list(main='Item $item_id',xlab="Test score",
                                                ylab="Probability",bty='l'),
                                   override=list(x = c(0,maxScore),y = 0:1, type="n"))

        plot.args$main = fstr(plot.args$main, list(item_id=i))
        plot.args$sub = fstr(plot.args$sub, list(item_id=i))

        do.call(plot, plot.args)

        draw_curtains(qnt)
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j], lwd=3)
        }
        prb = zR[ssI$first : ssI$last,]
        pte = bty(nrow(prb))
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j])
        }
      } # eol items
      
    } # eo not summate
    
  } # eo not overlay
}



print.rim <- function(x, ...){
  res = paste0('Parameters for the Rasch and Interaction Model\n\nitems: ',
               paste0(pull(x$inputs$ssI,"item_id"), collapse = ', '),
               '\n\n# use plot() for plotting the Rasch and Interaction Model or coef() for retreiving the parameters\n')
  cat(res)
  invisible(res)
}


# TODO: reparameterize pars for IM
coef.rim = function(object, ...) 
{
  x = object
  first = x$inputs$ssI$first
  last  = x$inputs$ssI$last
  dRM = x$est$bRM[-first]
  dRM = -log(dRM/dRM[1])
  dIM = x$est$bIM[-first]
  dIM = -log(dIM/dIM[1])
  OPCML_RM = toOPLM(x$inputs$ssIS$item_score, x$est$bRM, first, last)
  
  IS = tibble(item_id = x$inputs$ssIS$item_id[-first], item_score = x$inputs$ssIS$item_score[-first],
              beta_rasch = as.vector(OPCML_RM$delta), beta_IM = dIM)
  I = tibble(item_id = x$inputs$ssI$item_id, sigma = log(x$est$cIM), SE_sigma= x$est$se.c, fit_IM=x$est$fit.stats)
  
  inner_join(IS,I,by='item_id') %>% arrange(.data$item_id, .data$item_score) %>% as.data.frame()
}

