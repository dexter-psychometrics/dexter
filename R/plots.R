
##############################
my_layout <- function(npic, nr, nc) {
  if(npic==1) nr=nc=1
  nc = min(nc, 3)
  nc = min(nc, npic)
  nw = npic %/% nc + npic %% nc
  nr = min(nr, 3)
  nr = min(nr, nw)
  list(nr=nr, nc=nc)
}


# format string with named arguments
# @param f format string with named arguments prefixed by a dollar sign, formatting can be done with postfixing with : 
# @param arglist list with named arguments to interpolate in the format string. Use only alphanumerical characters in the names
# @return 
# formatted string
# @examples
# fstr('$bla, $b',list(bla='some string'))
# fstr('$bla:.1f, $b',list(bla=3.2423))
#
fstr = function(f, arglist)
{
  if(length(arglist)==0 || is.null(f)) return(f)
  f = gsub(paste0('\\$(?!(\\$|',paste0('(',names(arglist),')',collapse='|'),'))'),'$$', f, perl = TRUE)
  f = gsub('\\:(?!(\\.|\\:|(\\+\\.)))','::', f, perl=TRUE)
  rprintf(f,arglist)
}



###################################################
#' Distractor plot
#'
#' Produce a diagnostic distractor plot for an item
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, response, item_score
#' @param item The ID of the item to plot. A separate plot will be produced
#' for each booklet that contains the item, or an error message if the item ID
#' is not known. Each plot contains a non-parametric regression of each possible
#' response on the total score.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param ... further arguments to plot.
#' @details 
#' Customisation of title and subtitle can be done by using the arguments main and sub. 
#' These arguments can contain references to the variables item_id, booklet_id, item_position 
#' (only if dataSrc is a dexter db), pvalue, rit and rir. References are made by prefixing 
#' these variables with a dollar sign. Variable names can optionally be postfixed 
#' with a sprintf string, e.g.:
#'  
#' \code{distractor_plot(
#' db, main='item: $item_id', 
#' sub='Item-rest corr.: $rir:.2f')}
#' 
distractor_plot <- function(dataSrc, item, predicate = NULL, nc=1, nr=1, ...){  
  qtpredicate = eval(substitute(quote(predicate)))
  
  if(is.null(qtpredicate) & inherits(dataSrc,'DBIConnection'))
  {
    # cheat a little to make things faster
    qtpredicate = quote(booklet_id %in% booklets)
    env = new.env()
    env$booklets = dbGetQuery(dataSrc,'SELECT booklet_id FROM dxBooklet_design WHERE item_id=:item;',tibble(item=item))$booklet_id
  } else
  {
    env = caller_env()
  }
  
  # this uses the custom filter method on dx_resp_data objects
  respData = get_resp_data(dataSrc, qtpredicate = qtpredicate, extra_columns='response', env=env, summarised=FALSE) %>%
    filter(.data$item_id == item, .recompute_sumscores = FALSE )
  
  respData$x$response = factor(respData$x$response)

  if (nrow(respData$design) == 0) stop(paste("Item", item, "not found in dataSrc."))
  
  default.args = list(sub = "Pval: $pvalue:.2f, Rit: $rit:.3f, Rir: $rir:.3f", 
                      xlab = "Sum score", ylab = "Proportion", cex.sub = 0.8, xaxs="i")
  
  default.args$main = ifelse('item_position' %in% colnames(respData$design),
                             '$item_id: position $item_position in booklet $booklet_id',
                             '$item_id in booklet $booklet_id')
  
  user.args = list(...)

  rsp_counts = respData$x %>%
    group_by(.data$booklet_id, .data$response, .data$item_score, .data$sumScore) %>%
    summarise(n = n()) %>% 
    ungroup() 
  
  max_score = max(rsp_counts$item_score)
  # statistics by booklet
  stats = respData$x %>%
    group_by(.data$booklet_id) %>% 
    summarise(pvalue = mean(.data$item_score)/max_score, 
              rit = cor(.data$item_score, .data$sumScore), 
              rir = cor(.data$item_score, .data$sumScore - .data$item_score)) %>%
    ungroup() %>%
    inner_join(respData$design, by='booklet_id') %>%
    arrange(.data$booklet_id)

  npic = nrow(stats)
  ly = my_layout(npic, nr, nc)
  graphics::layout(matrix(1:(ly$nr * ly$nc), byrow = TRUE, 
                          ncol = ly$nc))
  
  labelz = levels(rsp_counts$response)
  
  rsp_counts = split(rsp_counts, rsp_counts$booklet_id)
  
  stats %>%
    rowwise() %>%
    do({
      st = as.list(.)
  
      y = rsp_counts[[as.character(st$booklet_id)]]
      
      plot.args = merge_arglists(user.args, 
                                  default=default.args,
                                  override=list(x = c(0,max(y$sumScore)), y = c(0,1), type="n"))
      
      plot.args$main = fstr(plot.args$main, st)
      plot.args$sub = fstr(plot.args$sub, st)

      do.call(graphics::plot, plot.args)

      bkl_scores = y %>% 
        group_by(.data$sumScore) %>% 
        summarise(n = sum(.data$n)) %>%
        ungroup() %>%
        arrange(.data$sumScore)
      
      dAll = density(bkl_scores$sumScore, n = 51, weights = bkl_scores$n/sum(bkl_scores$n))
      N = sum(bkl_scores$n)
      
      # this generates the legend and has the  side effect of drawing the lines in the plot
      lgnd = y %>% group_by(.data$response)  %>% do({
        dxi = density(.$sumScore, n = 51, weights = .$n/sum(.$n), 
                      bw = dAll$bw, from = min(dAll$x), to = max(dAll$x))
        yy = dxi$y/dAll$y * sum(.$n)/N
        k = match(.$response[1], labelz) + 1
        graphics::lines(dAll$x, yy, co = k, lw = 2)
        tibble(col = k, resp = paste0(.$response[1]," (", .$item_score[1], ")"))
      })
      graphics::legend("right", legend = as.character(lgnd$resp), 
                       lty = 1, col = lgnd$col, cex = 0.8, box.lty = 0)
    
      data.frame()
    })
  if(nc>1 || nr>1) par(mfrow=c(1,1))
  invisible(NULL)
}



#' Profile plot
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: 
#' person_id, item_id, item_score and the item_property and the covariate of interest.
#' @param item_property The name of the item property defining the domains. 
#' The item property should have exactly two distinct values in your data
#' @param covariate name of the person property/covariate used to create the groups. 
#' There will be one line for each distinct value.
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @param model "IM" (default) or "RM" where "IM" is the interaction model and 
#' "RM" the Rasch model. The interaction model is the default as it fits 
#' the data better or at least as good as the Rasch model.
#' @param x Which value of the item_property to draw on the x axis, if NULL, one is chosen automatically
#' @param ... further arguments to plot, many have useful defaults
#' @return Nothing interesting
#' @details 
#' Profile plots can be used to investigate whether two (or more) groups of respondents 
#' attain the same test score in the same way. The user must provide a  
#' (meaningfull) classification of the items in two non-overlapping subsets such that 
#' the test score is the sum of the scores on the subsets. 
#' The plot shows the probabilities to obtain 
#' any combinations of subset scores with thin gray lines indicating the combinations 
#' that give the same test score. The thick lines connect the most likely 
#' combination for each test score in each group.
#' When applied to educational test data, the plots can be used to detect differences in the 
#' relative difficulty of (sets of) items for respondents that belong to different 
#' groups and are matched on the test score. This provides a content-driven way to 
#' investigate differential item functioning. 
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#' covariates=list(gender="<unknown>"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' profile_plot(db, item_property='mode', covariate='gender')
#' 
#' close_project(db)
#' }
#' 
#' 
profile_plot <- function(dataSrc, item_property, covariate, predicate = NULL, model = "IM", x = NULL, ...) 
{
  if (model != "IM") model="RM"
  user.args = list(...)
  if(!inherits(dataSrc,'data.frame'))
  {
    item_property = tolower(item_property)
    covariate = tolower(covariate)
  }
  
  qtpredicate = eval(substitute(quote(predicate)))
  respData = get_resp_data(dataSrc, qtpredicate, extra_columns = covariate, 
                           extra_design_columns=item_property, env = caller_env()) 
  
  # make sure we have an intersection
  if(length(unique(respData$design$booklet_id)) > 1)
  {
    common_items = Reduce(intersect, split(respData$design$item_id, respData$design$booklet_id))
    if(length(common_items) == 0) stop('The intersection of the items across booklets is empty')
    
    respData$design = respData$design %>%
      distinct(.data$item_id, .keep_all=TRUE) %>%
      semi_join(tibble(item_id = common_items), by = 'item_id') %>%
      mutate(booklet_id='b')
    
    respData$x = respData$x %>%
      semi_join(respData$design, by='item_id') %>%
      group_by(.data$person_id) %>%
      mutate(sumScore = sum(.data$item_score), booklet_id = 'b') %>%
      ungroup()
  }
  
  if(nrow(respData$x) == 0) stop('no data to analyse')
  
  props = unique(respData$design[[item_property]])
  if(length(props) != 2)
    stop('this function needs an item_property with 2 unique values in your data')
  # order props to get the x and y axis right
  if(!is.null(x)) if(x %in% props) props = c(x, props[props!=x])
  
  
  # fit model
  models = by(respData$x, respData$x[[covariate]], function(rsp)
  {
    ssIS = rsp %>% group_by(.data$item_id, .data$item_score) %>%
      summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>%
      ungroup()
    
    ssIS = ssIS %>% 
      full_join(tibble(item_id = unique(ssIS$item_id), item_score = 0), by=c('item_id','item_score')) %>%
      arrange(.data$item_id, .data$item_score)
    
    ssIS[is.na(ssIS)] = 0
    
    ssI = ssIS %>% 
      group_by(.data$item_id) %>%
      summarise(nCat = n(), N = sum(.data$sufI), sufC = sum(.data$sufC), mx = max(.data$item_score)) %>%
      ungroup() %>%
      arrange(.data$item_id) %>% 
      mutate(first = cumsum(.data$nCat) - .data$nCat + 1, last = cumsum(.data$nCat))
    

    ssT = rsp %>% 
      group_by(.data$sumScore) %>% 
      summarise(N=n_distinct(.data$person_id)) %>%
      ungroup() %>%
      right_join(tibble(sumScore=0:sum(ssI$mx)), by='sumScore') %>%
      arrange(.data$sumScore)
    
    ssT[is.na(ssT)] = 0
    
    ss = list(il = ssI, sl = ssIS, tl = ssT)
    list(est=EstIM(ss), ss=ss)
  })
  
  
  tt = lapply(models, function(x)
  {
    A = c(1:nrow(x$ss$il))[x$ss$il$item_id %in% respData$design[respData$design[[item_property]]==props[1],]$item_id]
    B = c(1:nrow(x$ss$il))[x$ss$il$item_id %in% respData$design[respData$design[[item_property]]==props[2],]$item_id]
    SSTable(x, AB = list(A,B), model = model)
  })
  
  
  maxA = max(sapply(tt,function(x){ nrow(x$tbl)} ))-1
  maxB = max(sapply(tt,function(x){ ncol(x$tbl)} ))-1
  
  sg = data.frame(k=0:(maxA+maxB))
  
  default.args = list(asp=1, main="Profile plot", xlab=props[1], 
                      ylab=props[2],xlim=c(0,maxA),ylim=c(0,maxB))
  do.call(graphics::plot, 
          merge_arglists(user.args, 
                         default=default.args,
                         override=list(x=c(0,maxA), y=c(0,maxB),
                                       xaxs="i", type="n")))
  
  # The timolines
  k = maxA + maxB
  sg$y0 = pmin(maxB,sg$k)
  sg$x0 = sg$k - sg$y0
  sg$x1 = pmin(maxA,sg$k)
  sg$y1 = sg$k - sg$x1
  graphics::segments(sg$x0, sg$y0, sg$x1, sg$y1, col="gray")
  
  graphics::text(0:maxA,0,0:maxA,cex=.6,col="lightgray")
  graphics::text(maxA,1:maxB,(maxA+1:maxB),cex=.6,col="lightgray")
  
  
  for (i in seq_along(tt)) {
    ta = tt[[i]]$tbl
    y = tibble(
      value=as.vector(ta),
      Var1=as.integer(gl(nrow(ta),1,nrow(ta)*ncol(ta))),
      Var2=as.integer(gl(ncol(ta),nrow(ta),nrow(ta)*ncol(ta))),
      v = as.integer(gl(nrow(ta),1,nrow(ta)*ncol(ta))) + 
        as.integer(gl(ncol(ta),nrow(ta),nrow(ta)*ncol(ta)))
    )
    stp = y %>% 
      group_by(.data$v) %>%
      do(.[which.max(.$value),]-1)
    graphics::lines(stp$Var1, stp$Var2, col=i+1, lw=2)
  }
  
  graphics::legend("topleft", 
                   legend=names(tt), 
                   lty=1, col=1+(1:length(tt)),
                   cex=.7, 
                   box.lty=0)
  graphics::box()
}

