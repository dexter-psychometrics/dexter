
library(dplyr)

RcppArmadillo::armadillo_throttle_cores(1)


test_that('multiple b',{
  db = open_project(test_path('testdata/verbAggression.db'))
  
  ndr = 50
  f = fit_enorm(db,method='Bayes',nDraws=ndr)
  
  theta=seq(-3,3,by=.1)
  
  
  max_score = coef(f) |>
    group_by(item_id) |>
    summarise(m=max(item_score)) |>
    summarise(m=sum(m)) |>
    pull('m')
    
  
  pf = p_score(f, parms_draw='sample')(theta)
  
  expect_true(all(dim(pf) ==c(length(theta), max_score+1, ndr)), label='p_score dim correct')
  
  # all scores occur, so each row should have a single top and each column should have a single top 
  
  single_top = function(x)
  {
    w = which.max(x)
    !(is.unsorted(x[1:w]) || is.unsorted(x[length(x):w]))
  }
  expect_true(all(sapply(1:ndr, function(i)
    {
      all(apply(pf[,,i],1,single_top)) &&  all(apply(pf[,,i],2,single_top))
    })),label='pf single top over rows and columns')
  
  
  ii = information(f, parms_draw='sample')(theta)
  
  expect_true(all(apply(ii,2,single_top)), label='information topped')
  
  i=information(f, parms_draw='av')(theta)
  

  expect_gt(cor(i,rowMeans(ii),method='sp'), .99,
            label='information over average pars is approx equal to average information across pars')
  
  
  es =  expected_score(f, parms_draw='sample')(theta)
  

  expect_true(!any(apply(es,2,is.unsorted)),label='expected score increase over theta')
  
  
  x = r_score(f, parms_draw='sample')(theta)
  
  expect_true(all(dim(x) ==c(length(theta), nrow(get_items(f)), ndr)), label='dimensions multiple pars r_score') 
  
  skip_on_cran()
  
  # only local, chance elements
  expect_gt(mean(sapply(1:ndr,\(i) cor(rowSums(x[,,i]),theta))),.95,label='sim data high correlation sumscore and theta')

})


RcppArmadillo::armadillo_reset_cores()