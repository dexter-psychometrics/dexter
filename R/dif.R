

## produces a matrix of statistics for pairwise DIF
PairDIF_ <- function(par1,par2,cov1,cov2)
{
  labs=rownames(par1)
  D=kronecker(par2,t(par2),FUN="-")-kronecker(par1,t(par1),FUN="-") 
  var1=diag(cov1)
  var2=diag(cov2)
  S=(kronecker(var1,t(var1),FUN="+")-2*cov1)+(kronecker(var2,t(var2),FUN="+")-2*cov2)
  diag(S)=1
  D=D/sqrt(S)
  colnames(D)=labs; rownames(D)=labs
  return(D)
}

## produces a statistics for overall-DIF
OverallDIF_ <- function(par1,par2, cov1,cov2)
{
  r=1
  nI=length(par1)
  beta=par1-par2
  Sigma=cov1+cov2
  DIF_test=mahalanobis(beta[-r],rep(0,(nI-1)),Sigma[-r,-r])
  DIF_p=pchisq(DIF_test,(nI-1),lower.tail=FALSE)
  return(list(stat=DIF_test,df=nI-1, p=DIF_p))
}

###


#################################
bty = function (n, h = c(265, 75), c. = c(61, 66),
                l = c(25, 73), power = c(0.7, 1.742),
                fixup = TRUE, gamma = NULL, alpha = 1, ...)
{
  if (!is.null(gamma))
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L)
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- rep(c., length.out = 2L)
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, 0, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                       C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) *
                         rval), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}




#' Exploratory test for Differential Item Functioning
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param person_property Defines groups of persons to calculate DIF
#' @return An object of class \code{DIF_stats} holding statistics for
#' overall-DIF and a matrix of statistics for DIF in the relative position of
#' item-category parameters in the regular parameterization used e.g., by OPLM.
#' @details 
#' Tests for equality of relative item difficulties category locations across groups.
#' Supplements the confirmatory approach of the profile plot
#' 
#' @references 
#' Bechger, T. M. and Maris, G (2015); A Statistical Test for Differential Item Pair Functioning. 
#' Psychometrika. Vol. 80, no. 2, 317-340.
#' 
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", person_propertys=list(gender='unknown'))
#' add_booklet(db, verbAggrData, "agg")
#' dd = DIF(db,person_property="gender")
#' print(dd)
#' plot(dd)
#' 
#' close_project(db)
#' }
#' 
DIF = function(dataSrc, person_property, predicate=NULL) 
{

  qtpredicate = eval(substitute(quote(predicate)))
  
  if(!inherits(dataSrc,'data.frame'))
  {
    person_property = tolower(person_property)
  }
  
  columns = c('person_id','item_id','item_score', person_property)
  
  respData = get_responses_(dataSrc, qtpredicate, columns = columns, env = caller_env()) 
  
  if(nrow(respData) == 0) stop('no data to analyse')
  
  respData[[person_property]] = as.character(respData[[person_property]])
  
  # split into list
  respData = split(respData, respData[[person_property]])
  
  if(length(respData)!=2)
    stop('The person_property needs to have two unique values in your data to calculate DIF')
  
  common_scores = Reduce(dplyr::intersect, lapply(respData, function(x) x[,c('item_id','item_score')]))
  problems = Reduce(dplyr::union,
                    lapply(respData,
                           function(data) dplyr::setdiff(data[,c('item_id','item_score')],common_scores)))
  
  if(nrow(problems) > 0)
  {
    problems = tibble(item_id = unique(problems$item_id))
    warning(paste('the following items do not have the same score categories over both person_propertys and',
                  'have been removed from the analysis:', paste0(problems$item_id,collapse=', ')))
    respData = lapply(respData, function(x) x %>% anti_join(problems, by = 'item_id'))
  }

  ## 2. Estimate models with fit_enorm using CML
  models = lapply(respData, fit_enorm)
  #models = list(fit_enorm_(respData[[1]],env=as.environment(list())), 
  #              fit_enorm_(respData[[2]],env=as.environment(list())))
  


  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = OverallDIF_ (models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
                           models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)
  
  D = PairDIF_(models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
               models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)
  

  items = common_scores %>%
    filter(.data$item_score > 0) %>%
    arrange(.data$item_id, .data$item_score)

  ## 5. Report D and DIF_stats and inputs
  ou = list(DIF_overall = DIF_stats, DIF_pair = D, 
            group_labels = names(respData), items = items)
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
  yLabels <- rev(yLabels)
  DIF_pair <- DIF_pair[nrow(DIF_pair) : 1,]
  
  # Data Map
  oldpar = par(mar = c(6,8,2.5,2))
  on.exit({par(oldpar)},add=TRUE)
  do.call(image,
          merge_arglists(list(...),
                         override = list(x = 1:length(yLabels), y = 1:length(xLabels), z=t( DIF_pair),
                                         col=ColorRamp),
                         default = default.args))
  
  
  axis(1, at=1:length(yLabels), labels=yLabels, las= 3, cex.axis=0.6, hadj=1,padj=1)
  axis(2, at=1:length(xLabels), labels=xLabels, las=1, cex.axis=0.6, hadj=1,padj=1)
  
  #Color Scale
  par(mar=c(6,3,2,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  graphics::layout(1)

}