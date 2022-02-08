
####################################################
# Computes likelihood and test information for internal use
#
# For a vector of thetas it returns:
# l = a matrix (nbr of response cats * length of theta) of the likelihood or log-likelihood if log=TRUE
# I = a vector of the information function computed at each theta = sum(P'^2/PQ)
# J = a vector of the derivative of the information function at each theta
#
# The vector theta can be a set of quadrature points or estimates to compute their SE
#
# Note: can not deal with Inf or NA values in theta
IJ_ = function(b, a, first, last, theta, log=FALSE)
{
  nI = length(first)
  nT = length(theta)
  I = matrix(0, nT, nI)
  J = matrix(0, nT, nI)
  logFi = double(nT)

  a = as.integer(a)
  
  IJ_c(theta, b, a, as.integer(first-1L), as.integer(last-1L), I, J,logFi)

  scores = 0:sum(a[last])
  
  l = sweep(outer(scores,theta), 2, logFi, '-')
  if (!log) l = exp(l)
  return(list(I=rowSums(I), J=rowSums(J), l=l))
}


# simulate test-scores rather then response patterns. Adapted for inclusion zero
rscore = function(theta,b,a,first,last, cntr=NULL, use_b_matrix=FALSE)
{
  if(use_b_matrix)
  {
    if(is.null(cntr)) stop('use_b_matrix is true, need a counter')
    b = b[cntr(),]
  }
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  a = as.integer(a)
  
  sampleNRM2_test(theta, b, a, first, last)[,1,drop=TRUE]

}

rscore_item = function(theta,b,a,first,last)
{
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  a = as.integer(a)
  sampleNRM2_item(theta, b, a, first, last)
}



# Expected scores given one or more ability values theta
E_score = function(theta,b,a,first,last)
{
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  a = as.integer(a)
  
  Escore_C(theta, b, a, first, last)[,1,drop=TRUE]
}

# Estimate a single ability for a whole score distribution.
# testing for overdispersion
theta_score_distribution = function(b,a,first,last,scoretab)
{
  ms.a = sum(a[last])
  theta = 0
  np = sum(scoretab)
  escore = -1
  score = ((0:ms.a) %*% scoretab)[1,1,drop=TRUE]

  first = as.integer(first-1L)
  last = as.integer(last-1L)
  
  while (abs(escore-score)>1e-6)
  {
    escore = np * Escore_C(theta,b,a,first,last)
    theta = theta + log(score/escore)
  }
  return(theta)
}

# exp(score*tht) is liable to get infinite
# g is largish so only compounds the problem
# p has cols theta and rows scores


# Expected distribution given a vector theta
# return matrix, ncol=length(theta), nrow=nscores
# old
# pscore = function(theta, b, a, first, last)
# {
#   g = elsym(b, a, first, last)
#   score = 0:(length(g)-1)
#   
#   p = sapply(theta, function(tht) g * exp(score*tht))
#   sweep(p,2,colSums(p),`/`)
# }

pscore = function(theta, b, a, first, last)
{
  g = elsym(b, a, first, last)
  score = 0:(length(g)-1)
  p = sapply(theta, function(tht) log(g) + score*tht)
  
  exp(sweep(p,2,apply(p,2,logsumexp),`-`))
}


# vector of 0,1 indicating if a score is possible. element 1 is score 0
possible_scores = function(a, first, last)
  drop(possible_scores_C(as.integer(a), as.integer(first-1L), as.integer(last-1L)))


theta_MLE <- function(b,a,first,last, se=FALSE)
{
  a = as.integer(a)
  theta = theta_mle_sec(b, a, as.integer(first-1L), as.integer(last-1L))[,1,drop=TRUE]
  
  sem = NULL
  if (se)
  {
    # use r indexed first last for IJ
    f = IJ_(b,a,first,last, theta)
    sem = c(NA, 1/sqrt(f$I), NA)
  }

  return(list(theta = c(-Inf,theta,Inf), se=sem))
}

theta_WLE <- function(b,a,first,last, se=FALSE)
{
  a = as.integer(a)
  theta = theta_wle_sec(b, a, as.integer(first-1L), as.integer(last-1L))[,1,drop=TRUE]
  
  sem = NULL
  if (se)
  {
    # use r indexed first last for IJ
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  
  return(list(theta = theta, se=sem))
}




## EAP using Jeffrey's prior: aka propto sqrt(information)
# Uses a weighted average to integrate over a grid defined by:
# grid_from, grid_to and grid_length.
theta_jEAP = function(b, a, first,last, se=FALSE, grid_from=-6, grid_to=6, grid_length=101) 
{
  theta_grid = seq(grid_from, grid_to, length=grid_length)
  f = IJ_(b,a,first,last,theta_grid)
  prior=sqrt(f$I)
  w = sweep(f$l, 2, prior, '*')
  theta = apply(w, 1, function(x) weighted.mean(theta_grid, w=x))
  sem=rep(NA,length(theta))
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  return(list(theta=theta,se=sem))
}



# se is always returned (arg se is ignored)
theta_EAP_GH = function(b, a, first,last, se=TRUE, mu=0, sigma=4)
{
  nodes = quadpoints$nodes * sigma + mu
  weights = quadpoints$weights
  ps = t(pscore(nodes,b,a,first,last))
  theta_EAP_GH_c(ps,nodes,weights)
}


##### EAP based on npv (default 500) plausible values ####
# Uses recycling to get npv plausible values for each sum score.
# Options:
# Smoothing is an option to avoid non-monotonicity due to random fluctuations
# se is the option to calculate standard-errors or not
### In R-code:
# n=rep(npv,length(score))
# R=matrix(0,length(score),npv)
# while (any(n>0))
# {
#   atheta=rnorm(1,mu,sigma)
#   sm=rscore(b,a,atheta,first,last)+1
#   if (n[sm]>0)
#   {
#     R[sm,n[sm]]=atheta
#     n[sm]=n[sm]-1
#   }
# }
# The final argument allows EAPs based on A-weighted score, where A need not equal a.
####

  

theta_EAP = function(b, a, first, last, npv=500, mu=0, sigma=4, smooth=FALSE, se=FALSE, A=NULL)
{
  if(is.null(A)) A=a
  
  mx = sum(A[last])
  score = (0:mx)[as.logical(possible_scores(A,first,last))]
  tmp = pv_recycle(b,a,first,last,score,npv,mu,sigma,A=A)
  
  theta = rep(NA,(mx+1))
  theta[score+1]=rowMeans(tmp)
  
  # @Timo: ik denk niet dat smooth goed werkt als er NA scores zijn
  # misschien gewoon alleen een vector voor de bestaande scores gebruiken?
  if (smooth) 
  {
      score = 0:mx
      theta = predict(lm(theta ~ poly(score,7)))
  }
  sem=rep(NA,(mx+1))

  if (se) sem=apply(tmp,1,sd)
  return(list(theta=theta, se=sem))
}




## Wrapper to C function
# currently using mean_ElSym
ittotmat = function(b,c,a,first,last)
{
  ittotmat_C(b,as.integer(a),c,as.integer(first-1L), as.integer(last-1L))
}



## Wrapper to C function
# using Elsym
ittotmat0 = function(b,c,a,first,last,ps)
{
  ittotmat0_C(b,as.integer(a),c,as.integer(first-1L), as.integer(last-1L), as.integer(ps))
}


########################################################
## Score-by-score table. Currently using mean_ElSym
########################################################

# @param AB       list: two mutually exclusive subsets of items as indexes first/last/etc.
# @return         a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b
# @details        NA's indicate that a total scores was not possible given the weights
# if cIM is not null, the interaction model will be used
SSTable = function(b, a, first, last, AB, cIM=NULL)
{
  design = tibble(first = as.integer(first - 1L),
                  last = as.integer(last - 1L),
                  set = NA_character_)
  
  design$set[AB[[1]]] = 'A'
  design$set[AB[[2]]] = 'B'
  
  A = filter(design, .data$set=='A')
  B = filter(design, .data$set=='B')
  
  if(!is.null(cIM))
  {
    if(anyNA(design$set))
      stop('for IM the two subsets should make up the complete set')
    
    ic = cIM[rep(1:nrow(design), design$last - design$first + 1L)]
    ss_table_im_C(a, b, ic, design$first, design$last,
                  A$first, A$last, B$first, B$last)
    
  } else
  {
    design = filter(design, !is.na(.data$set))
    
    ss_table_enorm_C(a, b, design$first, design$last,
                             A$first, A$last, B$first, B$last)
  }
}

# Polynomial smoothing of the log-lambda's
smooth_log_lambda = function(log_lambda, degree, robust=TRUE)
{
  score_range = 0:(length(log_lambda)-1)
  degree = min(degree, sum(!is.na(log_lambda)))
  if (robust){
    qr = lmsreg(log_lambda ~ poly(score_range, degree, raw=TRUE))
  }else
  {
    qr = lm(log_lambda ~ poly(score_range, degree, raw=TRUE))
  }
  predict(qr, new=data.frame(score_range))
}


## Get the score distribution of a booklet from fit_enorm
#  based on a polynomial smoothing of the log-lambda's
#  Currently only implemented for CML
# TO DO: Implement for Bayes.
# Check e.g., plot(0:48,log(lambda),col="green"); lines(0:48,log_l_pr)
# coeficients beta = as.numeric(qr$coefficients)[-1]
# n.obs is the exact observed score distributions if CML
ENORM2ScoreDist <- function(b, a, lambda, first, last, degree=2) 
{
  
  log_l_pr = smooth_log_lambda(log(lambda), degree=degree)

  g = elsym(b,a,first,last)
  lambda[is.na(lambda)] = 0
  sc_obs = g*lambda
  sc_sm = g*exp(log_l_pr)
  
  data.frame(score    = 0:sum(a[last]),
             n.obs    = sc_obs, 
             n.smooth = sc_sm,
             p.obs    = sc_obs/sum(sc_obs),
             p.smooth = sc_sm/sum(sc_sm))
}






