

# colors derived from http://colorbrewer2.org
qcolors = function(n, user_colors=NULL)
{
  if(is.function(user_colors))
  {
    pal = user_colors(n)
  } else if(length(user_colors) > 1 || length(user_colors) == n)
  {
    pal = user_colors
  } else
  {
    pal = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999",
            "#BC80BD", "#CCEBC5","#FFED6F","#FB8072","#80B1D3","#FDB462","#8DD3C7","#FFFFB3","#BEBADA",
            "#B3DE69","#A50F15","#08306B","#00441B","#54278F","#fc4E2A","#525252","#66C2A4")
  }
  rep_len(pal,n)
}

dvcolors = function(n, user_colors=NULL)
{
  if(is.function(user_colors))
  {
    pal = user_colors(n)
  } else if(length(user_colors) > 1 || length(user_colors) == n)
  {
    pal = user_colors
  } else
  {
    pal = c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5',
            '#3288BD','#5E4FA2')
  }
  colorRampPalette(pal)(n)
}



lighten = function(color, factor = 0.5) {
  if ((factor > 1) || (factor < 0)) stop("factor needs to be within [0,1]")
  col = col2rgb(color)
  col = col + (255 - col)*factor
  col = rgb(t(col), maxColorValue=255)
  col
}


draw_curtains = function(qnt)
{
  if(!is.null(qnt))
  {
    usr = par('usr')
    rect(usr[1], usr[3], qnt[1], usr[2], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
    rect(qnt[2], usr[3], usr[2], usr[4], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
    # rect will sometimes cover the axis (redrawing the axis would probably solve this too)
    abline(h=usr[3])
    abline(v=usr[1])
  }
}


###################################################
#' Distractor plot
#'
#' Produce a diagnostic distractor plot for an item
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: person_id, item_id, response, item_score
#' and optionally booklet_id
#' @param item_id The ID of the item to plot. A separate plot will be produced
#' for each booklet that contains the item, or an error message if the item_id
#' is not known. Each plot contains a non-parametric regression of each possible
#' response on the total score.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param legend logical, whether to include the legend. default is TRUE
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param adjust factor to adjust the smoothing bandwidth respective to the default value
#' @param col vector of colors to use for plotting. The names of the vector can be responses. If the vector is not named, 
#' colors are assigned to the most frequent responses first.
#' @param ... further arguments to plot.
#' @return 
#' Silently, a data.frame of response categories and colors used. Potentially useful if you want to customize the legend or 
#' print it separately
#' @details 
#' Customization of title and subtitle can be done by using the arguments main and sub. 
#' These arguments can contain references to the variables item_id, booklet_id, item_position(if available),
#' pvalue, rit and rir. References are made by prefixing these variables with a dollar sign. Variable names may be postfixed 
#' with a sprintf style format string, e.g. 
#' \code{distractor_plot(db, main='item: $item_id', sub='Item rest correlation: $rir:.2f')}
#' 
distractor_plot = function(dataSrc, item_id, predicate=NULL, legend=TRUE, curtains=10, adjust=1, col=NULL, ...){  
  check_dataSrc(dataSrc)

  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  item_id = as.character(item_id)
  check_string(item_id)

  user.args = list(...); leg.args = list()
  if(length(names(user.args))>0)
  {
    leg.args = user.args[endsWith(names(user.args),'.legend')]
    names(leg.args) = gsub('\\.legend$','',names(leg.args))
    user.args = user.args[!endsWith(names(user.args),'.legend')]
  }
  
  iprop = list()
  if(is_db(dataSrc))
    iprop = as.list(dbGetQuery_param(dataSrc,'SELECT * FROM dxItems WHERE item_id= :item_id;', 
                                     tibble(item_id=item_id)))
  
  
  if(is.null(qtpredicate) && is_db(dataSrc))
  {
    # pre process a little to make things faster
  	booklets = dbGetQuery_param(dataSrc,
  	     'SELECT booklet_id FROM dxbooklet_design WHERE item_id=:item_id;', tibble(item_id=item_id)) |>
  		pull('booklet_id') |>
  	  sql_quote("'") |>
  		paste(collapse=',')
  	
    qtpredicate = sql(paste0("booklet_id IN(",booklets,")"), 'booklet_id')
  } 
  item = item_id
  respData = get_resp_data(dataSrc, qtpredicate = qtpredicate, extra_columns='response', env=env, summarised=FALSE) |>
    filter_rd(.data$item_id == item, .recompute_sumscores = FALSE )
  

  if (nrow(respData$design) == 0) 
    stop(paste("Item", item_id, "not found in dataSrc."))
  

  if('item_position' %in% colnames(respData$design))
  {
    ipos = filter(respData$design, .data$item_id==!!item_id)
  } else if(is_bkl_safe(dataSrc, qtpredicate, env) && is_db(dataSrc))
  {
    ipos = dbGetQuery(dataSrc, paste("SELECT booklet_id, item_position FROM dxbooklet_design WHERE item_id=", 
                                     sql_quote(item_id,"'")))
  } else if(inherits(dataSrc,'data.frame') && 'item_position' %in% colnames(dataSrc))
  {
    ipos = dataSrc |>
      filter(.data$item_id==!!item_id) |>
      distinct(.data$booklet_id,.keep_all=TRUE)
  } else
  {
    ipos = filter(respData$design, .data$item_id==!!item_id)
  }
  ipos=select(ipos,any_of(c('booklet_id','item_position')))
  
  if(inherits(dataSrc,'data.frame'))
     respData$x = mutate(respData$x,response=coalesce(as.character(.data$response),'<NA>'))

  default.args = list(sub = "Pval: $pvalue:.2f, Rit: $rit:.3f, Rir: $rir:.3f", 
                      xlab = "Sum score", ylab = "Proportion", cex.sub = 0.8, xaxs="i", bty="l")
  
  default.args$main = ifelse('item_position' %in% colnames(ipos),
                             '$item_id, pos. $item_position in booklet $booklet_id',
                             '$item_id in booklet $booklet_id')
  

  rsp_counts = respData$x |>
    count(.data$booklet_id, .data$response, .data$item_score, .data$booklet_score)

  max_score = max(rsp_counts$item_score)

  rsp_colors = rsp_counts |>
    group_by(.data$response) |>
    summarise(n = sum(.data$n)) |>
    ungroup() |>
    arrange(desc(.data$n), .data$response) 
  
  # cannot make this a list since a response can be an empty string
  if(is.null(names(col)))
  {
    rsp_colors$color = qcolors(nrow(rsp_colors),col)
  } else
  {
    rsp_colors$color = col[rsp_colors$response]
    rsp_colors$color[rsp_colors$response==''] = col[names(col)==''][1] # nog ff checken
    rsp_colors$color[is.na(rsp_colors$color)] = '#A9A9A9'
  }
  
  
  lapply(split(rsp_counts, as.character(rsp_counts$booklet_id)), function(y)
  {
    booklet = y$booklet_id[1]
    
    bkl_scores = y |> 
      group_by(.data$booklet_score) |> 
      summarise(n = sum(.data$n)) |>
      ungroup() |>
      arrange(.data$booklet_score)
    
    if(nrow(bkl_scores)>1)
    {
      N = sum(bkl_scores$n)
      
      labs = modifyList(iprop,
        list(pvalue = weighted.mean(y$item_score, y$n)/max_score,
                  rit = weighted_cor(y$item_score, y$booklet_score, y$n),
                  rir = weighted_cor(y$item_score, y$booklet_score - y$item_score, y$n),
                  n = N,
                  item_position = filter(ipos, .data$booklet_id==booklet)$item_position,
                  booklet_id=booklet,
                  item_id=item_id))
      
      plot.args = merge_arglists(user.args, 
                                 default=default.args,
                                 override=list(x = c(0,max(y$booklet_score)), y = c(0,1), type="n"))
      
      plot.args$main = fstr(plot.args$main, labs)
      plot.args$sub = fstr(plot.args$sub, labs)
      
      qua = curtains/200
      qnt = NULL
      if(qua>0 && qua<.5) {
        qnt = weighted_quantile(bkl_scores$booklet_score, bkl_scores$n, c(qua,1-qua))
      }
      
      do.call(plot, plot.args)
      draw_curtains(qnt)
      
      dAll = density(bkl_scores$booklet_score, n = 512, weights = bkl_scores$n/N, adjust=adjust,
                     from=min(bkl_scores$booklet_score),to=max(bkl_scores$booklet_score),warnWbw=FALSE)
      # has to be from=min(booklet_score), otherwise we can get division by 0
      lgnd = y |> 
        group_by(.data$response)  |> 
        do({
          k = rsp_colors[rsp_colors$response == .$response[1],]$color
          
          if(nrow(.)==1)
          {
            yval = .$n/sum(filter(y, .data$booklet_score==.$booklet_score)$n)
            points(.$booklet_score,yval,col=k,pch=16,xpd=TRUE)
          } else
          {
            dxi = density(.$booklet_score, weights = .$n/sum(.$n), n = 512,
                        bw = dAll$bw, from = min(dAll$x), to = max(dAll$x),warnWbw=FALSE)
            yy = dxi$y/dAll$y * sum(.$n)/N
            lines(dAll$x, yy, col = k, lw = 2)
          }
          tibble(col = k, label = paste0(.$response[1]," (", .$item_score[1], ")"))
        })
      
      if(legend && NROW(lgnd)>0)
      {
        do.call(graphics::legend,
                merge_arglists(leg.args, 
                               default=list(x="topleft", cex=.8, box.lty=0, bg='white',lwd=2,inset=c(0.01, 0)),
                               override=list(legend=lgnd$label,lty=1, col=lgnd$col)))
      }
    }
  })

  invisible(df_format(rsp_colors))
}


adjust_pp_scores = function(maxA, psbl,cex)
{
  w = max(strwidth(maxA,cex=cex), strheight('1',cex=cex)) + 0.05
  
  if(w>=1)
  {
    lag = psbl[1]
    for(i in 2:length(psbl))
    {
      if(psbl[i]-lag<w)
      {
        psbl[i] = NA_integer_
      } else
      {
        lag = psbl[i]
      }
    }
  }
  psbl[!is.na(psbl)]
}

pp_segments = function(maxA, maxB, psbl, col='lightgray',cex=0.6)
{
  psbl = adjust_pp_scores(maxA,psbl,cex)
  
  clip(0,maxA,0,maxB) 
  segments(0L, psbl, psbl, 0L, col=col, xpd=FALSE)
  
  h = strheight('1',cex=cex)
  
  y = par('usr')[3] + h
  psblx = psbl[psbl<=maxA]
  
  text(psblx, y, psblx, cex=cex,col=col, xpd=TRUE,pos=4,offset=0)
  
  psbly = psbl[psbl>maxA]
  text(maxA ,psbly-maxA+y, psbly,col=col,xpd=TRUE,cex=cex,pos=4,offset=0)
  
}


#' Profile plot
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: 
#' person_id, item_id, item_score and the item_property and the covariate of interest.
#' @param item_property The name of the item property defining the domains. 
#' The item property should have exactly two distinct values in your data
#' @param covariate name of the person property used to create the groups. 
#' There will be one line for each distinct value.
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @param model "IM" (default) or "RM" where "IM" is the interaction model and 
#' "RM" the Rasch model. The interaction model is the default as it fits 
#' the data better or at least as good as the Rasch model.
#' @param x Which category of the item_property to draw on the x axis, if NULL, one is chosen automatically
#' @param col vector of colors to use for plotting
#' @param col.diagonal color of the diagonal lines representing the testscores
#' @param ... further graphical arguments to plot. Graphical parameters for the legend can be postfixed with \code{.legend}
#' @details 
#' Profile plots can be used to investigate whether two (or more) groups of respondents 
#' attain the same test score in the same way. The user must provide a  
#' (meaningful) classification of the items in two non-overlapping subsets such that 
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
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules, ":memory:", 
#'                          person_properties=list(gender="unknown"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' profile_plot(db, item_property='mode', covariate='gender')
#' 
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
profile_plot = function(dataSrc, item_property, covariate, predicate = NULL, model = c("IM","RM"), x = NULL, 
                        col = NULL, col.diagonal='lightgray',...) 
{
  check_dataSrc(dataSrc)
  check_string(item_property)
  check_string(covariate)
  model = match.arg(model)
  
  user.args = list(...); leg.args = list()
  if(length(names(user.args))>0)
  {
    leg.args = user.args[endsWith(names(user.args),'.legend')]
    names(leg.args) = gsub('\\.legend$','',names(leg.args))
    user.args = user.args[!endsWith(names(user.args),'.legend')]
  }
  
  
  if(is_db(dataSrc))
  {
    item_property = dbValid_colnames(item_property)
    covariate = dbValid_colnames(covariate)
  }

  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  if(inherits(dataSrc,'matrix'))
	stop('profile_plot does not accept a matrix dataSrc')
  
  respData = get_resp_data(dataSrc, qtpredicate, extra_columns = covariate, env = env)  |>
	  intersection_rd()
  
  respData$x[[covariate]] = ffactor(respData$x[[covariate]])
  
    
  if(is_db(dataSrc) && item_property %in% dbListFields(dataSrc,'dxitems'))
  {
    domains = dbGetQuery(dataSrc, paste("SELECT item_id,", item_property, "FROM dxitems;")) |>
      semi_join(tibble(item_id = levels(respData$x$item_id)), by='item_id')
    
    domains$item_id = ffactor(domains$item_id, levels = levels(respData$x$item_id))
    
  } else
  {
    domains = distinct(dataSrc, .data$item_id, .data[[item_property]])
    
    if(nrow(domains) > n_distinct(domains$item_id))
      stop("some items belong to multiple domains, this is not allowed")
  }
  
  domains = arrange(domains, .data$item_id)
  
  props = unique(domains[[item_property]])
  if(length(props) != 2)
    stop('this function needs an item_property with 2 unique values in your data')
  
  # order props to get the x and y axis right
  if(!is.null(x) && x %in% props) 
    props = c(x, props[props!=x])
  
  setA = which(domains[[item_property]] == props[1])
  setB = which(domains[[item_property]] == props[2])             

  # fit interaction model and compute ssTable
  models = by_rd(respData, covariate, fit_inter_, regs=FALSE) 
  
  tt = lapply(models, function(m)
  {
      if(model=="IM")
      {
        SSTable(b=m$est$bIM, a=m$inputs$ssIS$item_score, 
                first=m$inputs$ssI$first, last=m$inputs$ssI$last, setA=setA, setB=setB, cIM_score = m$est$cIM_score)
      } else
      {
        SSTable(b=m$est$bRM, a=m$inputs$ssIS$item_score, 
                first=m$inputs$ssI$first, last=m$inputs$ssI$last, setA=setA, setB=setB)
      }
  })

  psbl = Reduce(union,lapply(models,function(m) m$est$possible_scores))
  
  maxA = max(sapply(tt, nrow)) - 1L
  maxB = max(sapply(tt, ncol)) - 1L

  default.args = list(main="Profile plot", xlab=props[1], 
                      ylab=props[2],xlim=c(0,maxA),ylim=c(0,maxB),bty='l')
  do.call(plot, 
          merge_arglists(user.args, 
                         default=default.args,
                         override=list(x=c(0,maxA), y=c(0,maxB),type="n")))
  
  seg.cex = 0.6 * ifelse('cex' %in% names(user.args), user.args$cex, 1)
  pp_segments(maxA,maxB,psbl,cex=seg.cex, col=col.diagonal)
  
  colors = qcolors(length(tt), col)
  
  for (i in seq_along(tt)) 
  {
    lns = sapply(psbl, function(s){
        indx = cbind(1L+(s:0),1L+(0:s))
        indx = indx[indx[,1]<=nrow(tt[[i]]) & indx[,2] <= ncol(tt[[i]]),, drop=FALSE]
        indx = indx[which.max(tt[[i]][indx]),, drop=FALSE]
        if(nrow(indx) > 0 && tt[[i]][indx] > 0)
          return(indx)
        rep(NA_integer_, 2)
      }) |>
      t() |>
      apply(2, '-', 1) 
    
    lns = lns[!(is.na(lns[,1]) | is.na(lns[,2])),]
    lines(lns,col=colors[i], lw=2)
  }

  do.call(legend,
          merge_arglists(leg.args, 
                         default=list(x="topleft", cex=.7, box.lty=0, bg='white',inset=0.01),
                         override=list(legend=names(tt),lty=1, col=colors)))

  invisible(NULL)
}


plot.prms = function(x, item_id=NULL, dataSrc=NULL, predicate=NULL, nbins=5, ci = .95, 
                      add=FALSE, col = 'black', col.model='grey80', ...)
{
  plot.enorm(x, item_id=item_id, dataSrc=dataSrc, predicate=predicate, nbins=nbins, ci = ci, 
             add=add, col = col, col.model=col.model, ...)
  
}

#' Plot for the extended nominal Response model
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_enorm
#' @param item_id which item to plot, if NULL, one plot for each item is made
#' @param dataSrc data source, see details
#' @param predicate an expression to subset data in dataSrc
#' @param nbins number of ability groups
#' @param ci confidence interval for the error bars, between 0 and 1. Use 0 to suppress the error bars.
#' Default = 0.95 for a 95\% confidence interval
#' @param add logical; if TRUE add to an already existing plot
#' @param col color for the observed score average
#' @param col.model color for the expected score based on the model
#' @param ... further arguments to plot
#' @return 
#' Silently, a data.frame with observed and expected values possibly useful to create a numerical fit measure.
#' @details
#' The standard plot shows the fit against the sample on which the parameters were fitted. If
#' dataSrc is provided, the fit is shown against the observed data in dataSrc. This may be useful 
#' for plotting the fit in different subgroups as a visual test for item level DIF. The confidence 
#' intervals denote the uncertainty about the predicted pvalues within the ability groups for the 
#' sample size in dataSrc (if not NULL) or the original data on which the model was fit.
#' 
#' @examples
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' db = start_new_project(verbAggrRules, ":memory:", 
#'   person_properties=list(gender=""))
#' 
#' add_booklet(db, verbAggrData, "agg")
#' 
#' f = fit_enorm(db)
#' 
#' plot(f, items="S1DoShout")
#' 
#' # side by side for two different groups
#' # (it is also possible to show two lines in the same plot 
#' # by specifying add=TRUE as an argument in the second plot)
#' 
#' par(mfrow=c(1,2))
#' 
#' plot(f,items="S1WantCurse",dataSrc=db, predicate = gender=='Male', 
#'   main='men - $item_id')
#' 
#' plot(f,items="S1WantCurse",dataSrc=db, predicate = gender=='Female', 
#'   main='women - $item_id')
#' 
#' close_project(db)
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
#' @method plot enorm
#' 
plot.enorm = function(x, item_id=NULL, dataSrc=NULL, predicate=NULL, nbins=5, ci = .95, 
                     add=FALSE, col = 'black', col.model='grey80', ...)
{
  check_num(nbins,'integer',.length=1, .min=2)
  dots = list(...)
  
  if(inherits(x,"mst_enorm"))
  {
    m = x$inputs$method
    x$inputs = x$mst_inputs
    x$inputs$method = m
  }
  
  if(is.null(item_id))
  {
    if('items' %in% names(dots))
    {
      # common typo
      item_id = dots$items
      dots$items = NULL
    } else
    {
      item_id = x$inputs$ssI$item_id
    }
  }
  if(length(setdiff(item_id,x$inputs$ssI$item_id))>0)
  {
    message('The following items were not found in your fit object')
    print(setdiff(item_id,x$inputs$ssI$item_id))
    stop('unknown item',call.=FALSE)
  }
  
  if(!is.null(dataSrc))
  {
    check_dataSrc(dataSrc)
    qtpredicate = eval(substitute(quote(predicate)))
    env = caller_env()
    respData = get_resp_data(dataSrc, qtpredicate, env=env, retain_person_id=FALSE,
                             parms_check=filter(x$inputs$ssIS, .data$item_id %in% local(item_id)))
    
    if(length(setdiff(as.character(item_id), levels(respData$design$item_id)))>0)
    {
      message('The following items were not found in dataSrc')
      print(setdiff(as.character(item_id), levels(x$design$item_id)))
      stop('unknown item',call.=FALSE)
    }
    
    x$abl_tables = list(mle = filter(ability_tables(x,design=respData$design, method='MLE',parms_draw='average'),is.finite(.data$theta)) )
    if(!is.factor(x$abl_tables$mle$booklet_id))
      x$abl_tables$mle$booklet_id = factor(x$abl_tables$mle$booklet_id,levels=levels(respData$design$booklet_id))
    
    x$inputs$plt = get_sufStats_nrm(respData, check_sanity=FALSE)$plt
  }
  
  
  #many plots
  if(length(item_id) > 1)
  {
    out = lapply(item_id, function(itm) do.call(plot, append(list(x=x, item_id=itm, nbins=nbins, ci=ci), dots)))
    names(out) = as.character(item_id)
    
    return(invisible(out))
  }
  # for dplyr
  item_id_ = as.character(item_id)
  
  expf = expected_score(x, items = item_id)
  
  max_score = x$inputs$ssIS |>
    filter(.data$item_id == item_id_) |>
    pull(.data$item_score) |>
    max()
  
  plt = x$inputs$plt |>
    filter(.data$item_id==item_id_) |>
    inner_join(x$abl_tables$mle, by=c('booklet_id','booklet_score')) |>
    mutate(abgroup = weighted_ntile(.data$theta, .data$N, nbins = nbins)) |>
    group_by(.data$abgroup) |>
    summarize(gr_theta = weighted.mean(.data$theta,.data$N), avg_score = weighted.mean(.data$meanScore,.data$N), n=sum(.data$N)) |>
    ungroup() |>
    mutate(expected_score = expf(.data$gr_theta))
  
  rng = max(plt$gr_theta) - min(plt$gr_theta)
  rng = c(min(plt$gr_theta)-.5*rng/nbins,
          max(plt$gr_theta)+.5*rng/nbins)
  
  plot.args = merge_arglists(dots,
                             default=list(bty='l',xlab = expression(theta), ylab='score',main=item_id,
                                          lwd=par('lwd')),
                             override=list(x = rng,y = c(0,max_score), type="n"))
  
  plot.args$main = fstr(plot.args$main, list(item_id=item_id))
  plot.args$sub = fstr(plot.args$sub, list(item_id=item_id))
  
  if(!add) do.call(plot, plot.args)
  
  plot(expf,from = rng[1], to=rng[2], col=col.model, add=TRUE,lwd=plot.args$lwd)
  
  plt$outlier = FALSE
  
  if(!is.null(ci) && !is.na(ci) && ci !=0)
  {
    if(ci>1 && ci<100)
      ci = ci/100
    
    if(ci<0 || ci >= 1)
      stop('confidence interval must be between 0 and 1')
    
    qnt = abs(qnorm((1-ci)/2))
    
    I=information(x, items = item_id)
    plt = plt |>
      mutate(se = sqrt(I(.data$gr_theta)/.data$n),
             conf_min = pmax(.data$expected_score - qnt*.data$se,0),
             conf_max = pmin(.data$expected_score + qnt*.data$se,max_score)) |>
      mutate(outlier = .data$avg_score < .data$conf_min | .data$avg_score > .data$conf_max)
    
    suppressWarnings({
      arrows(plt$gr_theta, plt$conf_min, 
             plt$gr_theta, plt$conf_max, 
             lwd=plot.args$lwd,
             length=0.05, angle=90, code=3, col=col.model)})
  } 
  
  lines(plt$gr_theta,plt$avg_score,col=col,lwd=plot.args$lwd)  
  points(plt$gr_theta, plt$avg_score, 
         bg = if_else(plt$outlier, qcolors(1), coalesce(plot.args$bg,'transparent')), 
         pch = coalesce(dots$pch,21), col=col)
  invisible(df_format(plt))
}




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
#' @method plot inter
#'
plot.inter = function(x, items=NULL, summate=TRUE, overlay=FALSE,
                    curtains=10, show.observed=TRUE, ...){
  all_items = as.character(x$inputs$ssI$item_id)
  if(!is.null(items))
  {
    if(length(setdiff(items, all_items)) > 0) 
      stop_(format_plural('Items[s] not found: %s',setdiff(items, all_items)))
  } else
  {
    items = all_items
  }
  user.args = list(...)
  
  qua = curtains/200
  qnt=NULL
  if(qua>0 && qua<.5) {
    qnt = weighted_quantile(x$inputs$scoretab$booklet_score, x$inputs$scoretab$N, c(qua,1-qua))
  }
  
  npic = length(items)
  
  max_score = ncol(x$est$itrRM)-1L
  scores = x$est$possible_scores
  zRM = x$est$itrRM[,scores+1L]
  zIM = x$est$itrIM[,scores+1L]
  
  
  if (overlay && !summate) overlay=FALSE
  if (overlay) 
  {
    # only summate possible
    # do the Rasch model
    #
    z = zRM[row.names(zRM) %in% items,]
    items = row.names(z)
    maxy = max(z[,ncol(z)])
    
    plot.args = merge_arglists(user.args,
                               default=list(main="$model model",xlab="Test score",
                                            ylab="Item score", bty='l'),
                               override=list(x = c(0,max_score),y= c(0,maxy), type="n"))
    
    plot.args$main = fstr(plot.args$main, list(model='Rasch'))
    plot.args$sub = fstr(plot.args$sub, list(model='Rasch'))
    
    do.call(plot, plot.args)
    draw_curtains(qnt)
    
    for (i in 1:npic) 
      lines(scores,z[i,]) # the actual lines
    lx = sample(scores, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      points(lx[i], z[i,lx[i]+1], col="white", cex=1.6, pch=19)
      text(lx[i], z[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    # do the Interaction model
    #
    z = zIM[row.names(zIM) %in% items,]
    
    plot.args = merge_arglists(user.args,
                               default=list(main="$model model",xlab="Test score",
                                            ylab="Item score", bty='l'),
                               override=list(x = c(0,max_score),y= c(0,maxy), type="n"))
    
    plot.args$main = fstr(plot.args$main, list(model='Interaction'))
    plot.args$sub = fstr(plot.args$sub, list(model='Interaction'))
    
    do.call(plot, plot.args)
    draw_curtains(qnt)
    
    for (i in 1:npic) 
      lines(scores, z[i,])
    lx = sample(scores, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      points(lx[i], z[i,lx[i]+1], col="white", cex=1.6, pch=19)
      text(lx[i], z[i,lx[i]+1], items[i], col=1, cex=.6)
    }
    
    # end of overlay
    
  } else {
    # not overlay: do many plots
    
    if (summate) {
      # for each item in turn, do both models (with summation), and plot
      
      for (i in items) {
        maxy = zRM[row.names(zRM)==i,ncol(zRM)]
        
        plot.args = merge_arglists(user.args,
                                   default=list(main='Item $item_id',xlab="Test score",
                                                ylab="Item score", bty='l'),
                                   override=list(x = c(0,max_score),y = c(0,maxy), type="n"))
        
        plot.args$main = fstr(plot.args$main, list(item_id=i))
        plot.args$sub = fstr(plot.args$sub, list(item_id=i))
        
        do.call(plot, plot.args)
        
        draw_curtains(qnt)
        if(show.observed) {
          plt = filter(x$inputs$plt, .data$item_id == i)
          points(plt$booklet_score, plt$meanScore, col="coral",pch=20)
        }
        lines(scores, zIM[row.names(zIM)==i,], col="gray80", lwd=3)
        lines(scores, zRM[row.names(zRM)==i,])
      }
      
    } else {
      zI = x$est$ctrIM[,scores+1L]
      zR = x$est$ctrRM[,scores+1L]
      # for each item in turn, similar but possibly multiline and coloured
      for (i in items) {
        ssI = filter(x$inputs$ssI, .data$item_id==i)
        prb = zI[ssI$first : ssI$last, ]
        
        
        pte = lighten(qcolors(nrow(prb)),.6)
        
        plot.args = merge_arglists(user.args,
                                   default=list(main='Item $item_id',xlab="Test score",
                                                ylab="Probability",bty='l'),
                                   override=list(x = c(0,max_score),y = 0:1, type="n"))
        
        plot.args$main = fstr(plot.args$main, list(item_id=i))
        plot.args$sub = fstr(plot.args$sub, list(item_id=i))
        
        do.call(plot, plot.args)
        
        draw_curtains(qnt)
        for (j in 1:nrow(prb)) {
          lines(scores, prb[j,], col=pte[j], lwd=3)
        }
        prb = zR[ssI$first : ssI$last,]
        pte = qcolors(nrow(prb))
        for (j in 1:nrow(prb)) {
          lines(scores, prb[j,], col=pte[j])
        }
      } # eo items
      
    } # eo not summate
    
  } # eo not overlay
}
