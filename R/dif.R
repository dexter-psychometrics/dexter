

## produces a matrix of statistics for pairwise DIF
PairDIF_ <- function(beta1, beta2, acov.beta1, acov.beta2)
{
  labs = rownames(beta1)
  DR = kronecker(beta2,t(beta2),FUN="-")-kronecker(beta1,t(beta1),FUN="-") 
  var1 = diag(acov.beta1)
  var2 = diag(acov.beta2)
  S = (kronecker(var1,t(var1),FUN="+")-2*acov.beta1)+(kronecker(var2,t(var2),FUN="+")-2*acov.beta2)
  diag(S) = 1
  D = DR/sqrt(S)
  colnames(D) = labs; rownames(D)=labs
  colnames(DR) = labs; rownames(DR)=labs
  return(list(D=D, Delta_R=DR))
}

## produces a statistics for overall-DIF
# beta1 and beta1 are both mean centered and do not contain the zero category. 
# In general, this works if both sets of parameters have the same normalization.
# If a reference item is set to zero, r should be its index.
OverallDIF_ <- function(beta1, beta2, acov1, acov2)
{
  r = 1
  nI = length(beta1)
  d_beta = beta1-beta2
  Sigma = acov1+acov2
  DIF_test = mahalanobis(d_beta[-r],rep(0,(nI-1)),Sigma[-r,-r])
  DIF_p = pchisq(DIF_test,(nI-1),lower.tail=FALSE)
  return(list(stat=DIF_test,df=nI-1, p=DIF_p))
}


#' Exploratory test for Differential Item Functioning
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param person_property Defines groups of persons to calculate DIF
#' @return An object of class \code{DIF_stats} holding statistics for
#' overall-DIF and a matrix of statistics for DIF in the relative position of
#' item-category parameters in the beta-parameterization where they represent 
#' locations on the ability scale where adjacent categories are equally likely. 
#' If there is DIF, the function `plot` can be used to produce an image of the pairwise DIF statistics.
#' @details 
#' Tests for equality of relative item/category difficulties across groups.
#' Supplements the confirmatory approach of the profile plot.
#' 
#' @references 
#' Bechger, T. M. and Maris, G (2015); A Statistical Test for Differential Item Pair Functioning. 
#' Psychometrika. Vol. 80, no. 2, 317-340.
#' 
#' 
#' @seealso A plot of the result is produced by the function \code{\link{plot.DIF_stats}}
#' 
#' @examples
#' db = start_new_project(verbAggrRules, ":memory:", person_properties=list(gender='unknown'))
#' add_booklet(db, verbAggrData, "agg")
#' dd = DIF(db,person_property="gender")
#' print(dd)
#' plot(dd)
#' str(dd)
#' 
#' close_project(db)
#' 
DIF = function(dataSrc, person_property, predicate=NULL) 
{
  
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  check_dataSrc(dataSrc)
  check_string(person_property)
  
  if(is_db(dataSrc))
    person_property = tolower(person_property)
  
  pb = get_prog_bar(retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  ## 1. get data
  respData = get_resp_data(dataSrc, qtpredicate, extra_columns = person_property, env = env) 
  
  respData$x[[person_property]] = ffactor(as.character(respData$x[[person_property]]))
  
  if(nlevels(respData$x[[person_property]]) != 2)
    stop('The person_property needs to have two unique values in your data to calculate DIF')
  
  
  ## 2. Estimate models with fit_enorm using CML
  models = by(respData, person_property, fit_enorm)
  
  ## 3. get the intersection
  common_items = models[[1]]$inputs$ssIS %>%
    full_join(models[[2]]$inputs$ssIS, by=c('item_id','item_score') ,suffix=c('_1','_2')) %>%
    group_by(.data$item_id) %>%
    filter(!(anyNA(.data$sufI_1) | anyNA(.data$sufI_2))) %>%
    ungroup() %>%
    distinct(.data$item_id)
  
  if(nrow(common_items) != nrow(models[[1]]$inputs$ssI) || nrow(common_items) != nrow(models[[2]]$inputs$ssI))
  {
    cat('\n')
    message('Some items were excluded because they do not appear in both datasets and/or do not have the same score categories in both datasets.')
    models = lapply(models, function(m){
      ii = m$inputs$ssIS %>%
        filter(.data$item_score > 0L) %>%
        mutate(i = row_number()) %>%
        semi_join(common_items, by='item_id') %>%
        pull(.data$i)
        
      m$est$beta = drop(m$est$beta)[ii]
      m$est$acov.beta = m$est$acov.beta[ii,ii]
      return(m)
    })
  }
  
  
  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = OverallDIF_ (models[[1]]$est$beta, models[[2]]$est$beta, 
                           models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  DIF_mats = PairDIF_(models[[1]]$est$beta, models[[2]]$est$beta, 
                      models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  items = models[[1]]$inputs$ssIS %>%
    semi_join(common_items,by='item_id') %>%
    filter(.data$item_score > 0) %>%
    select(.data$item_id, .data$item_score) %>%
    arrange(.data$item_id, .data$item_score) %>%
    mutate(item_id=as.character(.data$item_id))
  
  
  ## 5. Report D and DIF_stats and inputs
  ou = list(DIF_overall = DIF_stats, DIF_pair = DIF_mats$D, Delta_R = DIF_mats$Delta_R, 
            group_labels = names(models), items = items)
  class(ou) = append('DIF_stats', class(ou))
  return(ou)
}




print.DIF_stats <- function(x, ...)
{
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  tmp = specify_decimal(x$DIF_overall$p,3)
  if (tmp=="0.000") tmp="< 0.0006"
  cat(paste0("Test for DIF:"," Chi-square = ", as.character(round(x$DIF_overall$stat, digits=3)),
             ", df = ", 
             as.character(x$DIF_overall$df),
             ", p = ", tmp))  
}


#' plot method for pairwise DIF statistics
#' 
#' @param x object produced by \code{\link{DIF}}
#' @param items character vector of item id's for a subset of the plot. Useful if you have many items. 
#' If NULL all items are plotted.
#' @param itemsX character vector of item id's for the X axis
#' @param itemsY character vector of item id's for the Y axis
#' @param alpha significance level used to color the plot (two sided)
#' @param ... further arguments to plot
#' 
#' @references
#' Feskens, R., Fox, J. P., & Zwitser, R. (2019). Differential item functioning in PISA due to mode effects. 
#' In Theoretical and Practical Advances in Computer-based Educational Measurement (pp. 231-247). Springer, Cham.
#' 
#' @details
#' Plotting produces an image of the matrix of pairwise DIF statistics. 
#' The statistics are standard normal deviates and colored to distinguish significant from non-significant values.
#' If there is no DIF, a proportion alpha off the cells will be colored significant by chance alone.
#'      
plot.DIF_stats = function(x, items = NULL, itemsX = items, itemsY = items, alpha =.05,...)
{
  if(is.null(itemsX)) itemsX = sort(unique(x$items$item_id))
  if(is.null(itemsY)) itemsY = sort(unique(x$items$item_id))
  
  if(length(setdiff(c(itemsX, itemsY), x$items$item_id)) > 0)
  {
    cat('items not found in DIF object:\n')
    print(setdiff(c(itemsX, itemsY), x$items))
    stop('some of the item_ids you specified are not present in the DIF object')
  }
  
  x$items = x$items %>%
    mutate(rn = row_number())
  
  itemsX = x$items %>%
    inner_join(tibble(item_id = itemsX, ord = 1:length(itemsX)), by='item_id') %>%
    arrange(.data$ord)
  
  itemsY = x$items %>%
    inner_join(tibble(item_id = itemsY, ord = 1:length(itemsY)), by='item_id') %>%
    arrange(.data$ord)
  
  DIF_pair = abs(x$DIF_pair[itemsX$rn, itemsY$rn])
  
  
  if(nrow(distinct(x$items,.data$item_id)) == nrow(x$items))
  {
    yLabels = pull(itemsY, 'item_id')
    xLabels = pull(itemsX, 'item_id')
  } else
  {
    yLabels = paste(itemsY$item_id, itemsY$item_score)
    xLabels = paste(itemsX$item_id, itemsX$item_score)
    
  }
  
  max_ = max(x$DIF_pair)
  default.args = list(main = paste(x$group_labels[1],'vs.',x$group_labels[2]),
                      axes=FALSE, zlim=c(0,max_),xlab='',ylab='',useRaster=TRUE)
  
  oldpar = par(no.readonly = TRUE)
  on.exit({par(oldpar)},add=TRUE)
  
  graphics::layout(matrix(c(1,1,2,0),2,2), widths=c(7,1))
  
  qn = qnorm(1-alpha/2)
  
  breaks = seq(0, qn, .05)
  col = colorRampPalette(c('white','lightblue'))(length(breaks)-1)
  if(max_ > qn)
  {
    breaks_b = seq(qn, min(20,max_), .05)
    breaks_b[length(breaks_b)] = max_ + 21
    col = c(col,colorRampPalette(c('gold1','red2'))(length(breaks_b)))
    breaks = c(breaks, breaks_b)
  }
  
  # Reverse Y axis
  # yLabels <- rev(yLabels)
  xLabels <- rev(xLabels)
  DIF_pair <- DIF_pair[nrow(DIF_pair) : 1,]
  
  # Data Map
  par(mar = c(6,8,2.5,2))
  
  
  do.call(image,
          merge_arglists(list(...),
                         override = list(x = 1:length(yLabels), y = 1:length(xLabels), z=t( DIF_pair),
                                         col=col,breaks=breaks),
                         default = default.args))
  
  
  axis(1, at=1:length(yLabels), labels=yLabels, las=3, cex.axis=0.6, hadj=1,padj=0.5)
  axis(2, at=1:length(xLabels), labels=xLabels, las=1, cex.axis=0.6, hadj=1,padj=0.5)
  
  #Color Scale
  try({
    par(mar=c(2,2,2.5,2))
    image(1, seq(0, max(qn,min(20,max_)), by=.05),
          matrix(seq(0, max(qn,min(20,max_)), by=.05),nrow=1),
          col=col,
          breaks=breaks,
          xlab="",ylab="",
          xaxt="n", axes=FALSE)
    axis(2, at=0:min(20,max_),lwd=0,lwd.ticks=1,las=2,cex.axis=0.8)
    
  }, silent=TRUE)
  
  graphics::layout(1)
  invisible(NULL)
}




DIF_heatmap = function(x, items = NULL, itemsX = items, itemsY = items, alpha =.05,...)
{
  if(is.null(itemsX)) itemsX = sort(unique(x$items$item_id))
  if(is.null(itemsY)) itemsY = sort(unique(x$items$item_id))
  
  if(length(setdiff(c(itemsX, itemsY), x$items$item_id)) > 0)
  {
    cat('items not found in DIF object:\n')
    print(setdiff(c(itemsX, itemsY), x$items))
    stop('some of the item_ids you specified are not present in the DIF object')
  }
  
  x$items = x$items %>%
    mutate(rn = row_number())
  
  itemsX = x$items %>%
    inner_join(tibble(item_id = itemsX, ord = 1:length(itemsX)), by='item_id') %>%
    arrange(.data$ord)
  
  itemsY = x$items %>%
    inner_join(tibble(item_id = itemsY, ord = 1:length(itemsY)), by='item_id') %>%
    arrange(.data$ord)
  
  DIF_pair = abs(x$DIF_pair[itemsX$rn, itemsY$rn])
  
  
  if(nrow(distinct(x$items,.data$item_id)) == nrow(x$items))
  {
    yLabels = pull(itemsY, 'item_id')
    xLabels = pull(itemsX, 'item_id')
  } else
  {
    yLabels = paste(itemsY$item_id, itemsY$item_score)
    xLabels = paste(itemsX$item_id, itemsX$item_score)
  }
  
  dimnames(DIF_pair) = list(xLabels, yLabels)
  max_ = max(x$DIF_pair)
  qn = qnorm(1-alpha/2)
  breaks = seq(0, qn, length=50)
  col = colorRampPalette(c('white','lightblue'))(length(breaks)-1)
  if(max_ > qn)
  {
    breaks_b = seq(qn+0.01, max_, length=50)
    col = c(col,colorRampPalette(c('gold1','red2'))(length(breaks_b)))
    breaks = c(breaks, breaks_b)
  }
  
  pheatmap::pheatmap(DIF_pair, colors=col, breaks=breaks, ...)  
  
  
  invisible(NULL)
}
