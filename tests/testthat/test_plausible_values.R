context('test plausible_values')

library(dplyr)


RcppArmadillo::armadillo_throttle_cores(1)

test_that('populations work',{
  db = open_project('../verbAggression.db')
  
  #artificially create two overlapping booklets
  f2 = fit_enorm(db, (gender=='Male' & item_position<=16) | 
                     (gender=='Female' & item_position>8))
  
  expect_true(n_distinct(f2$inputs$design$booklet_id)==2)
  
  set.seed(123)
  pv=plausible_values(db, f2, covariates='gender')
  
  expect_true(mean(pv[pv$gender=='Male',]$PV1) > mean(pv[pv$gender=='Female',]$PV1))
  
  
  
  abl = ability(db,f2,method='WLE')
  
  test = inner_join(abl,pv,by='person_id')
  
  expect_gt(cor(test$PV1,test$theta), .9,'verb agg correlation ability and plausible value should be larger than .9')
  
  dbDisconnect(db)
})
  
RcppArmadillo::armadillo_reset_cores()