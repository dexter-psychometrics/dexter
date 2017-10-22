

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
#' @param covariate Person covariate name for subsetting
#' @return An object of class \code{DIF_stats} holding statistics for
#' overall-DIF and a matrix of statistics for DIF in the relative position of
#' item-category parameters in the regular parameterization used e.g., by OPLM.
#' @details 
#' Tests for DIF as described in Bechger and Maris (2014; A Statistical Test 
#' for Differential Item Pair Functioning. Psychometrica). Supplements the 
#' confirmatory approach of the profile plot
#'
#' @examples
#' \dontrun{
#' db = start_new_project(verbAggrRules, "verbAggression.db", 
#' covariates=list(gender='unknown'))
#' add_booklet(db, verbAggrData, "agg")
#' dd = DIF(db,covariate="gender")
#' print(dd)
#' plot(dd)
#' 
#' close_project(db)
#' }
#' 
DIF = function(dataSrc, covariate, predicate=NULL) 
{
  ## 1. Interpret input.. much like beginning of profile plot
  # Check whether there are 2 groups 
  # Check whether predicate (if not NULL) does not leave us without data
  if(!inherits(dataSrc,'data.frame'))
  {
    covariate = tolower(covariate)
  }
  
  columns = c('person_id','item_id','item_score', covariate)
  
  qtpredicate=eval(substitute(quote(predicate)))
  respData = get_responses_(dataSrc, qtpredicate, columns = columns, env = caller_env()) 
  if(nrow(respData) == 0) stop('no data to analyse')
  
  # split into list
  respData = by(respData, respData[[covariate]], identity)
  
  if(length(respData)!=2)
    stop('The covariate needs to have two unique values in your data to calculate DIF')
  
  common_scores = Reduce(dplyr::intersect, lapply(respData, function(x) x[,c('item_id','item_score')]))
  problems = Reduce(dplyr::union,
                    lapply(respData,
                           function(data) dplyr::setdiff(data[,c('item_id','item_score')],common_scores)))
  
  if(nrow(problems) > 0)
  {
    problems = tibble(item_id = unique(problems$item_id))
    warning(paste('the following items do not have the same score categories over both covariates and',
                  'have been removed from the analysis:', paste0(problems$item_id,collapse=', ')))
    respData = lapply(respData, function(x) x %>% anti_join(problems, by = 'item_id'))
  }
  
  ## 2. Estimate models with fit_enorm using CML
  models = lapply(respData, fit_enorm)
  
  ## 3. Make sure parameters pertain to same items-responses in the same order
  # This should always be correct, since fit_enorm orders on items and scores and
  # I made sure in the prelim that both models have the same items and score catregories
  
  ## 4. Call overallDIF_ and PairDIF_
  DIF_stats = OverallDIF_ (models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
                           models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)
  
  D = PairDIF_(models[[1]]$est$beta.cml, models[[2]]$est$beta.cml, 
               models[[1]]$est$acov.cml, models[[2]]$est$acov.cml)
  
  
  ## 5. Report D and DIF_stats and inputs
  ou = list(DIF_overall = DIF_stats, DIF_pair = D, 
            groups = names(respData), items = unique(common_scores$item_id))
  class(ou) = append('DIF_stats', class(ou))
  return(ou)
}


print.DIF_stats <- function(x, ...)
{
  specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
  tmp = specify_decimal(x$DIF_overall$p,3)
  if (tmp=="0.000") tmp="< 0.0006"
  print(paste0("Test for DIF:"," Chi-square = ", as.character(round(x$DIF_overall$stat, digits=3)),
               ", df = ", 
               as.character(x$DIF_overall$df),
               ", p = ", tmp))  
}

plot.DIF_stats = function(x, ...)
{
  D = x #had to make the parameter name x instead of D for CRAN rules
  x=D$DIF_pair
  yLabels = rownames(x)
  xLabels = colnames(x)
  min_=min(x)
  max_=max(x)
  default.args = list(main = paste(D$groups, collapse = ' vs. '),
                      axes=FALSE, zlim=c(min_,max_),xlab='',ylab='')
  
  graphics::layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  
  tmp = rainbow(256)[1:128]
  ColorRamp=c(tmp, tmp[128:1])
  ColorLevels <- seq(min(x), max(x), length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  oldpar = par(mar = c(6,8,2.5,2))
  do.call(image,
          merge_arglists(list(...),
                         override = list(x = 1:length(xLabels), y = 1:length(yLabels), z=t(x),
                                         col=ColorRamp),
                         default = default.args))
  
  
  axis(1, at=1:length(xLabels), labels=xLabels, las= 3, cex.axis=0.6)
  axis(2, at=1:length(yLabels), labels=yLabels, las=1,
       cex.axis=0.6)
  
  #Color Scale
  par(mar=c(6,3,2,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  graphics::layout(1)
  par(oldpar)
  
}