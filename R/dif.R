

## produces a matrix of statistics for pairwise DIF
PairDIF_ <- function(beta1,beta2,acov.beta1,acov.beta2)
{
  labs=rownames(beta1)
  DR=kronecker(beta2,t(beta2),FUN="-")-kronecker(beta1,t(beta1),FUN="-") 
  var1=diag(acov.beta1)
  var2=diag(acov.beta2)
  S=(kronecker(var1,t(var1),FUN="+")-2*acov.beta1)+(kronecker(var2,t(var2),FUN="+")-2*acov.beta2)
  diag(S)=1
  D=DR/sqrt(S)
  colnames(D)=labs; rownames(D)=labs
  colnames(DR)=labs; rownames(DR)=labs
  return(list(D=D, Delta_R=DR))
}

## produces a statistics for overall-DIF
OverallDIF_ <- function(beta1,beta2, acov1,acov2)
{
  r=1
  nI=length(beta1)
  d_beta=beta1-beta2
  Sigma=acov1+acov2
  DIF_test=mahalanobis(d_beta[-r],rep(0,(nI-1)),Sigma[-r,-r])
  DIF_p=pchisq(DIF_test,(nI-1),lower.tail=FALSE)
  return(list(stat=DIF_test,df=nI-1, p=DIF_p))
}

# to do: not refer to oplm in help
# to do: coef for dif_stats?

#' Exploratory test for Differential Item Functioning
#'
#'
#' @param dataSrc Data source: a connection to a dexter database or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param person_property Defines groups of persons to calculate DIF
#' @return An object of class \code{DIF_stats} holding statistics for
#' overall-DIF and a matrix of statistics for DIF in the relative position of
#' item-category parameters in the regular parameterization used e.g., by OPLM.
#' @details 
#' Tests for equality of relative item/category difficulties across groups.
#' Supplements the confirmatory approach of the profile plot
#' 
#' @references 
#' Bechger, T. M. and Maris, G (2015); A Statistical Test for Differential Item Pair Functioning. 
#' Psychometrika. Vol. 80, no. 2, 317-340.
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
  
  if(!inherits(dataSrc,'data.frame'))
    person_property = tolower(person_property)

  respData = get_resp_data(dataSrc, qtpredicate, extra_columns = person_property, env = env) 

  respData$x[[person_property]] = ffactor(as.character(respData$x[[person_property]]))
  
  if(nlevels(respData$x[[person_property]]) != 2)
    stop('The person_property needs to have two unique values in your data to calculate DIF')
  
  item_scores = respData$x %>%
    distinct(.data[[person_property]], .data$item_id, .data$item_score) %>%
    count(.data$item_id, .data$item_score)
  
  if(!all(item_scores$n==2))
  {
    message('The following items do not have the same score categories over both person_properties and\n',
            'have been removed from the analysis:')
    
    remove = item_scores %>%
      filter(.data$n<2) %>%
      distinct(.data$item_id)
    
    print(remove$item_id)
    
    respData = respData %>%
      anti_join(remove, by='item_id', .recompute_sumscores = TRUE)
  }
  
  ## 2. Estimate models with fit_enorm using CML
  models = by(respData, person_property, fit_enorm)
  
  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = OverallDIF_ (models[[1]]$est$beta, models[[2]]$est$beta, 
                           models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  DIF_mats = PairDIF_(models[[1]]$est$beta, models[[2]]$est$beta, 
                      models[[1]]$est$acov.beta, models[[2]]$est$acov.beta)
  
  items = models[[1]]$inputs$ssIS %>%
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



#' plot method for DIF
#' 
#' @param x object produced by DIF
#' @param items character vector of item id's for a subset of the plot. Useful if you have many items. 
#' If NULL all items are plotted.
#' @param itemsX character vector of item id's for the X axis
#' @param itemsY character vector of item id's for the Y axis
#' @param ... further arguments to plot
#'    
plot.DIF_stats = function(x, items = NULL, itemsX = items, itemsY = items, ...)
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
  
  DIF_pair = x$DIF_pair[itemsX$rn, itemsY$rn]
  
  
  if(nrow(distinct(x$items,.data$item_id)) == nrow(x$items))
  {
    yLabels = pull(itemsY, 'item_id')
    xLabels = pull(itemsX, 'item_id')
  } else
  {
    yLabels = paste(itemsY$item_id, itemsY$item_score)
    xLabels = paste(itemsX$item_id, itemsX$item_score)
    
  }
  
  min_ = min(x$DIF_pair) # keep color range from complete object
  max_ = max(x$DIF_pair)
  default.args = list(main = paste(x$group_labels[1],'vs.',x$group_labels[2]),
                      axes=FALSE, zlim=c(min_,max_),xlab='',ylab='')
  
  graphics::layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  
  tmp = rainbow(256)[1:128]
  ColorRamp=c(tmp, tmp[128:1])
  ColorLevels <- seq(min(x$DIF_pair), max(x$DIF_pair), length=length(ColorRamp))
  
  # Reverse Y axis
  # yLabels <- rev(yLabels)
   xLabels <- rev(xLabels)
   DIF_pair <- DIF_pair[nrow(DIF_pair) : 1,]
  
  # Data Map
  oldpar = par(mar = c(6,8,2.5,2))
  on.exit({par(oldpar)},add=TRUE)
  do.call(image,
          merge_arglists(list(...),
                         override = list(x = 1:length(yLabels), y = 1:length(xLabels), z=t( DIF_pair),
                                         col=ColorRamp),
                         default = default.args))
  
  
  axis(1, at=1:length(yLabels), labels=yLabels, las= 3, cex.axis=0.6, hadj=1,padj=0.5)
  axis(2, at=1:length(xLabels), labels=xLabels, las=1, cex.axis=0.6, hadj=1,padj=0.5)
  
  #Color Scale
  par(mar=c(6,3,2,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  graphics::layout(1)

}