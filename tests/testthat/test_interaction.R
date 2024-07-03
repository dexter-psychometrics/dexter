context('check interaction model')

library(dplyr)

RcppArmadillo::armadillo_throttle_cores(1)


test_that('interaction model parameters are stable over simulation',{
  
  set.seed(123)
  db = open_project(test_path('verbAggression.db'))
  f = fit_inter(db)
  ts = get_testscores(db)
  close_project(db)
  
  simdat = r_score_IM(f, rep(ts$booklet_score,10))
  
  expect_true(all(rowSums(simdat) == rep(ts$booklet_score,10)),
              label='IM simulation returns sumscores')
  
  
  g = fit_inter(simdat)  
  
  f = coef(f)  
  g = coef(g)

  expect_gt(cor(f$beta_IM,g$beta_IM), 0.95, label='IM sim beta correlates >.95 true beta')
  
  i = seq(1,nrow(f),2)
  
  expect_gt(cor(f$sigma[i],g$sigma[i]), 0.9, label='IM sim sigma correlates >.9 true sigma')
    
})


RcppArmadillo::armadillo_reset_cores()