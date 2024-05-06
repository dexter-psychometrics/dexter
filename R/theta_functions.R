
#' Functions of theta
#' 
#' returns information function, expected score function, score simulation function, or score distribution 
#' for a single item, an arbitrary group of items or all items
#' 
#' @param parms object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and, 
#' depending on parametrization, a column named either beta/delta, eta or b
#' @param items vector of one or more item_id's. If NULL and booklet_id is also NULL, all items in parms are used
#' @param booklet_id id of a single booklet (e.g. the test information function), if items is not NULL this is ignored
#' @param parms_draw when the item parameters are estimated with method "Bayes" (see: \code{\link{fit_enorm}}), 
#' parms_draw specifies whether to use a sample (a different item parameter draw for each output column) or the posterior mean
#' of the item draws. Alternatively, it can be an integer specifying a specific draw. It is ignored when parms is not estimated Bayesianly.
#' @param by_chain add a separate dimension to the output denoting the chains in the Bayesian simulation. Useful for advanced users 
#' to check autocorrelation and mixing on substantive outcome measures. Defaults to FALSE. Ignored if parms is not estimated Bayesianly or parms_draw is an integer. 
#' @return Each function returns a new function which accepts a vector of theta's. These return the following values: 
#' \describe{
#' \item{information}{an equal length vector with the information estimate at each value of theta.}
#' \item{expected_score}{an equal length vector with the expected score at each value of theta}
#' \item{r_score}{a matrix with length(theta) rows and one column for each item containing simulated scores based on theta. 
#' To obtain test scores, use rowSums on this matrix}
#' \item{p_score}{a matrix with length(theta) rows and one column for each possible sumscore containing the probability of 
#' the score given theta}
#' }
#' 
#' @examples
#' 
#' \dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}
#' 
#' db = start_new_project(verbAggrRules,':memory:')
#' add_booklet(db,verbAggrData, "agg")
#' p = fit_enorm(db)
#' 
#' # plot information function for single item
#' 
#' ifun = information(p, "S1DoScold")
#' 
#' plot(ifun,from=-4,to=4)
#' 
#' # compare test information function to the population ability distribution
#' 
#' ifun = information(p, booklet="agg")
#' 
#' pv = plausible_values(db,p)
#' 
#' op = par(no.readonly=TRUE)
#' par(mar = c(5,4,2,4))
#' 
#' plot(ifun,from=-4,to=4, xlab='theta', ylab='test information')
#' 
#' par(new=TRUE)
#' 
#' plot(density(pv$PV1), col='green', axes=FALSE, xlab=NA, ylab=NA, main=NA)
#' axis(side=4)
#' mtext(side = 4, line = 2.5, 'population density (green)')
#' 
#' par(op)
#' close_project(db)
#' 
#' \dontshow{ RcppArmadillo::armadillo_reset_cores()}
#' 
information = function(parms, items=NULL, booklet_id=NULL, parms_draw = c('average','sample'), by_chain=FALSE)
{
  theta_function(parms, items=items, booklet=booklet_id, parms_draw=parms_draw, what='information', by_chain=by_chain)
}

#' @rdname information
expected_score = function(parms, items=NULL, booklet_id=NULL, parms_draw = c('average','sample'),by_chain=FALSE)
{
  theta_function(parms, items=items, booklet=booklet_id, parms_draw=parms_draw, what='expected', by_chain=by_chain)
}

#' @rdname information
r_score = function(parms, items=NULL, booklet_id=NULL, parms_draw = c('average','sample'),by_chain=by_chain)
{
  theta_function(parms, items=items, booklet=booklet_id, parms_draw=parms_draw, what='sim', by_chain=FALSE)
}

#' @rdname information
p_score = function(parms, items=NULL, booklet_id=NULL, parms_draw = c('average','sample'))
{
  theta_function(parms, items=items, booklet=booklet_id, parms_draw=parms_draw, what='pmf')
}



theta_function = function(parms, items=NULL, booklet=NULL, parms_draw=c('average','sample'), 
                          what=c('information','expected','sim','pmf'),
                          by_chain=FALSE)
{
  what = match.arg(what)
  if(is.factor(items)) items = as.character(items)
  check_character(items,nullable=TRUE)
  check_string(booklet,name='booklet_id',nullable=TRUE)

  # TO DO: chains and parms_draw only fully implemented for expected yet
  
  if(is.numeric(parms_draw)) parms_draw = as.integer(parms_draw)
  else parms_draw = match.arg(parms_draw)
  
  parms = simplify_parms(parms, draw=parms_draw, by_chain=by_chain)
  
  if(!is.null(items))
  {
    fl = parms$items |> 
      filter(.data$item_id %in% items)
    
    if(nrow(fl) != length(items))
    {
      message('unknown items:')
      print(sort(setdiff(items,fl$item_id)))
      stop('Some items were not found, see output')
    }
  } else if(!is.null(booklet))
  {
    fl = parms$design |> 
      filter(.data$booklet_id %in% booklet) |>
      distinct(.data$item_id, .keep_all=TRUE) |>
      arrange(.data$first)
    
    if(length(intersect(booklet,fl$booklet_id))<length(booklet))
    {
      stop('unknown booklet')
    }
  } else fl = parms$items
  a = parms$a
  b = parms$b
  multiple_b = !is.null(dim(b)) && ncol(b)>1
  if(parms_draw=='average') by_chain = FALSE
  if(by_chain && multiple_b)
  {
    chain_index = parms$chain_index
    mn = min(sapply(chain_index,length))
    chain_index = lapply(chain_index,function(i) i[1:mn])
  }
  
  rm(parms)

  #output
  
  if(what=='information')
  {
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta)))
        stop('theta may not contain nan/NA values')  
      
      w = is.finite(theta)
      if(multiple_b)
      {
        res = matrix(0,nrow = length(theta), ncol=nrow(b))
        for(i in 1:row(b))
          res[w,i] = IJ_(b[i,],a,fl$first, fl$last, theta[is.finite(theta)])$I
      } else
      {
        res = double(length(theta))
        res[w] = IJ_(b,a,fl$first, fl$last, theta[is.finite(theta)])$I
      }
      
      # extremely large values overflow to NaN, recover as 0
      res[is.nan(res)] = 0
      res
    }
    class(out) = append('inf_func',class(out))
    
  } else if(what == 'expected')
  {
    max_score = sum(a[fl$last])
    out = function(theta)
    {
      check_num(theta)
      if(any(is.na(theta) | is.nan(theta))) 
        stop('theta may not contain nan/NA values') 
      
      w = is.finite(theta)
      if(multiple_b && by_chain)
      {
        res = array(0,dim=c(length(theta), length(chain_index[[1]]), length(chain_index)))
        for(i in seq_along(chain_index))
        {
          for(j in 1:mn)
          {
            #cat(sprintf("%i:%i\n",i,indx))
            res[w,j,i] = E_score(theta[w],  
                                    b=b[chain_index[[i]][j],], a=a, 
                                    first=fl$first, last=fl$last)
          }
        }    
        #cat('done\n')
        res[!w & theta > 0,,] = max_score
      } else if(multiple_b)
      {
        res = matrix(0,nrow = length(theta), ncol=nrow(b))
        for(i in 1:nrow(b))
          res[w,i] = E_score(theta[w],  
                           b=b[i,], a=a, 
                           first=fl$first, last=fl$last)
        res[!w & theta > 0,] = max_score
      } else
      {
        res = double(length(theta))
      
        res[w] = E_score(theta[w],b=b, a=a, first=fl$first, last=fl$last)
        res[!w & theta > 0] = max_score
      }
      
      
      # extremely large values of theta overflow to NaN (small values undeflow to zero, which is fine)
      res[is.nan(res)] = max_score
      res
    }
    class(out) = append('exp_func',class(out))
  } else if(what=='sim')
  {
    # TO DO
    if(multiple_b) stop('not implementd')
    out = function(theta)
    {
      res = rscore_item(theta,b=b,a=a,first = fl$first, last = fl$last)
      colnames(res) = fl$item_id
      res
    }
    class(out) = append('sim_func',class(out))
  } else if(what=='pmf')
  {
    # TO DO
    if(multiple_b) stop('not implementd')
    out = function(theta)
    {
      res = pscore(theta,b=b,a=a,first = fl$first, last = fl$last)
      t(res)
    }
    class(out) = append('pmf_func',class(out))
  }
  
  out
}

print.inf_func = function(x,...) cat('Information function: I(theta)\n')
print.exp_func = function(x,...) cat('Conditional expected score function: E(X_i|theta)\n')
print.sim_func = function(x,...) cat('function to simulate item scores: (x_i1, ..., x_ip) ~ ENORM(theta)\n')
print.pmf_func = function(x,...) cat('Conditional score distribution function: P(x_+|theta)\n')



